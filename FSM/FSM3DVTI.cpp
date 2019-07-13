// Fast sweeping eikonal solver for 3D VTI medium
// Original C++ code by Tariq Alhalifah
// Mex wrapper by Ivan Abakumov
// Freie Universität Berlin
// e-mail: abakumov_ivan@mail.ru
// 20th of March 2018
// version: 0.1
//
//////////////////////////////////////////////////////
//c++  libraries

//c++  libraries
#include "mex.h"
#include <math.h>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <matrix.h>
#include <string.h>
#include <complex.h> 
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <numeric>
#include <iterator>
#include <vector>
#include <random>
#include <iostream>
#include <chrono>
#include <stdio.h>
#include <omp.h>

using namespace std;


///////////////////////////////////////////////////////////////////////////
////////////////////    GLOBAL VARIABLES     //////////////////////////////
///////////////////////////////////////////////////////////////////////////

#define info 0

#define PI 3.14159265358979323846

#define SF_HUGE 1e20

#include <rsf.h>
#include "fastvti.h"

#include "fastvti.h"


#define AT(x,y,z) ((z)+n3*(y)+n23*(x))
/*#define AT(z,y,x) ((z)+n1*(y)+n12*(x))*/





class GridClass{
    public:
        int    nx;
        int    ny;
        int    nz;
        int    nt;
        double x0;
        double dx;
        double y0;
        double dy;
        double z0;
        double dz;
        double t0;
        double dt;
        double mx;
        double my;
        double mz;
        double mt;
        double xx;
        double yy;
        double zz;
        double tt;
};

///////////////////////////////////////////////////////////////////////////
///////////////          FUNCTIONS          /////////////////////////
///////////////////////////////////////////////////////////////////////////


static double tp;
static float a, dd[8], rd1,rd2,rd3,d1,d2,d3;
static int nm, n23, n, n1, n2, n3, jj2, jj3;




static void update2 (int p1, int p2, int p3, float* tj, unsigned char* mj, float s, float rsv, float eta)
{
    float b, c, t1=0., t2=0., u, d3r=d3*s*rsv;
    float d3rr,c1,b1,g,gg,ggg,gggg,ggggg,ff;
    float A0,A1,A2,A3,A4,aaa1,bbb1,m1,m2,m3,S1,S2;
    double tp1,tp2,tp3,den=0.,tp4,tpA,tpB;
    unsigned int k, i;
    double aaa,bbb,ccc,ddd,ddd1;
    unsigned int jjj;

    b = c = 0.; i = k = 0;
    ddd=1.; jjj = 0;
    aaa = bbb = ccc = 0.; bbb1 = 0.; aaa1 = 0.;
    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.;
    tp1 = 0.; tp2 = 0.; tp3 = 0.; tp4 = 0.;

    if ((p1 > 0   ) && *(mj-n23) && ((t1 = *(tj-n23)) < *tj)) i |= 0x01;
    if ((p1 < n1-1) && *(mj+n23) && ((t2 = *(tj+n23)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; } 
    if (i & 0x01) {
	ddd=d1;
	if ((jjj & 0x01) && (p1 < n1-2) && *(mj+2*n23)) { 
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj+2*n23); t1 /= 3.;
	} else if ((p1 > 1) && *(mj-2*n23)) {   
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj-2*n23); t1 /= 3.;
	}
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1;
	i ^= 0x01; k |= 0x01;
    }
    jjj =0;

    if ((p2 > 0   ) && *(mj-n3) && ((t1 = *(tj-n3)) < *tj)) i |= 0x01;
    if ((p2 < n2-1) && *(mj+n3) && ((t2 = *(tj+n3)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; } 
    if (i & 0x01) {
	ddd=d2;
	if ((jjj & 0x01) && (p2 < n2-2) && *(mj+2*n3)) { 
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj+2*n3); t1 /= 3.; 
	} else if ((p2 > 1) && *(mj-2*n3)) {   
	    ddd *= 2.25; t1 *= 4.; t1 -= *(tj-2*n3); t1 /= 3.;
	}
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1;
	i ^= 0x01; k |= 0x02;
    }
    jjj =0;
    if(p2==jj2)
	sf_warning("p1=%d p2=%d p3=%d n1=%d n2=%d n3=%d",p1,p2,p3,n1,n2,n3);

    if(k){
	tp = (bbb +sqrt (SF_MAX(bbb*bbb-aaa*(ccc-s),0.)))/aaa;
	den = tp*aaa-bbb;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	tp1 = -(ccc-tp*(2.*bbb-tp*aaa))*den;
	tp2 = (tp1*(2.*(bbb-tp*aaa)-0.5*tp1*aaa))*den;
	tp3 = -(tp1*tp1*aaa-tp2*(2.*(bbb-tp*aaa)-tp1*aaa))*den;
	tp4 = -(tp2*aaa*(2*tp1+0.5*tp2)-tp3*(2.*(bbb-tp*aaa)-tp1*aaa))*den;
    }
    if(p2==jj2)
	sf_warning("tp=%f tp1=%f tp2=%f v=%f vv=%f eta=%f t1=%f den=%f",tp,tp1,tp2,1/sqrt(s),sqrt(rsv),eta,t1,den);

    ddd1=d3;
    if ((p3 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) i |= 0x01;
    if ((p3 < n3-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2; jjj ^= 0x01; }  
    if (i & 0x01) {
	if ((i & 0x02) && (!(i & 0x01) || t2 > t1)) t1 = t2; 
	ddd=d3r;
	if ((jjj & 0x01) && (p3 < n3-2) && *(mj+2)) { 
	    ddd *= 2.25; ddd1 *= 2.25; t1 *= 4.; t1 -= *(tj+2); t1 /= 3.; 
	} else if ((p3 > 1) && *(mj-2)) {   
	    ddd *= 2.25; ddd1 *= 2.25; t1 *= 4.; t1 -= *(tj-2); t1 /= 3.;
	}
	d3rr = -ddd1*ccc*rsv;
	u = d3rr*t1; b1 = bbb+u; c1 = ccc+u*t1;
	g = bbb+t1*aaa; gg = aaa; ggg = bbb*t1*t1;
	gggg = aaa*t1*t1; ggggg = t1*bbb; ff = aaa+d3rr;
	aaa += ddd; u = ddd*t1; bbb += u; ccc += u*t1; bbb1 = bbb; aaa1 = aaa;
	i ^= 0x01; k |= 0x04;

	ccc += -s;
	tp = (bbb +sqrt (SF_MAX(bbb*bbb-aaa*ccc,0.)))/aaa;
	den = tp*aaa1-bbb1;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	m1 = rsv*ddd1*(gggg+4.*ggggg);
	m2 = ff-m1+tp*ddd1*rsv*(3.*g-2.*tp*gg);
	m3 = b1-tp*m2-rsv*ddd1*ggg;
	tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(-(c1-2.*tp*b1+tp*tp*ff+tp*tp*tp*ddd1*rsv*
					  (2.*g-tp*gg)-tp*tp*m1+rsv*ddd1*tp*2.*ggg)*den),SF_ABS(tp1));
	tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS((2.*tp1*m3-0.5*tp1*tp1*aaa1)*den),SF_ABS(tp2));
	tp3 = SF_SIG(tp3)*SF_MIN(SF_ABS(-(tp1*tp1*(ff-m1+6.*rsv*ddd1*tp*
						   (g-tp*gg))-2.*tp2*m3+tp1*tp2*aaa1)*den),SF_ABS(tp3));
	tp4 = SF_SIG(tp4)*SF_MIN(SF_ABS(-(2.*tp1*(tp1*tp1*rsv*ddd1*(g-2.*tp*gg)+tp2*(ff-m1+6.*rsv*tp*ddd1*(g-tp*gg)))+
					  tp2*tp2*0.5*aaa1-2.*tp3*m3+tp1*tp3*aaa1)*den),SF_ABS(tp4));
    }
    jjj =0;

    if (!k) return;
    
    A0=tp; A1=A0+tp1*eta; A2=A1+tp2*eta*eta; A3=A2+tp3*eta*eta*eta; A4=A3+tp4*eta*eta*eta*eta;
    if(p2==jj2)
	sf_warning("tp=%f tp1=%f tp2=%f A1=%f A2=%f t1=%f den=%f",tp,tp1,tp2,tp4,A1,A2,t1,den);
    if(SF_ABS(tp1)>0.00000001 || SF_ABS(tp2)>0.00000001){
	den = tp1-tp2*eta;
	den = (SF_ABS(den)<0.000000000001 ? SF_SIG(den)*1000000000000. : 1./den);
	tpA = (tp*tp1-(tp*tp2-tp1*tp1)*eta)*den;
    } else
	tpA = tp;

  
    if(SF_ABS(tp1)>0.00000001 || SF_ABS(tp2)>0.00000001 || SF_ABS(tp3)>0.00000001 || SF_ABS(tp4)>0.00000001){
	den = A1+A3-2*A2;
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	S1 = (A3*A1-A2*A2)*den;
	den = A2+A4-2*A3;
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	S2 = (A4*A2-A3*A3)*den;
	den = tpA+S2-2*S1;
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	tpB = (S2*tpA-S1*S1)*den;
    } else
	tpB=tp;
    tp = tpA;
    if(p2==jj2)
	sf_warning("tp=%f den=%f k=%d j1=%d j2=%d j3=%d",tp,den,k,p1,p2,p3);
    if(tp<0 || tp >10){
	sf_warning("tp=%f tp1=%f den=%f k=%d j1=%d j2=%d j3=%d",tp,tp1,den,k,p1,p2,p3);
	exit(0);
    }
 

    if (tp < *tj) {
	*tj = tp;
	if (*mj == FMM_OUT) {
	    nm--; 
	    *mj = FMM_FRONT; 
	    sf_pqueue_insert (tj);
	} 
    }
}















///////////////////////////////////////////////////////////////////////////

/* Numerically solving the differential source linear equation */
void fastvti (float* time                /* time */, 
	          float* slow                   /* NMO velocity squared */, 
	          float* slowv                   /* vertical velocity squared */, 
	          float* eta                   /*  eta */, 
	          int* in                    /* in/front/out flag */, 
	          bool* plane                /* if plane source */,
	          int   nn3,  int nn2,  int nn1 /* dimensions */,
	          float o3,float o2,float o1 /* origin */,
	          float dd3,float dd2,float dd1 /* sampling */,
	          float ss3,float ss2,float ss1 /* source */,
	          int   b3,  int b2,  int b1 /* box around the source */,
	          int order                  /* accuracy order (1,2,3) */)
{
    float xs[3], *pt, *ttime;
    float sl = 0., svl = 0., etal = 0., sh;
    float s2, sv2, sint2, cost2, dist2, ss;
    int i, i1, i2, i3, j1, j2, j3, nh, j;
    unsigned char *mask, *pm;
    
    xs[2] = (ss1-o1)/dd1; d3 = dd1;  n3 = nn1;
    xs[1] = (ss2-o2)/dd2; d2 = dd2;  n2 = nn2;
    xs[0] = (ss3-o3)/dd3; d1 = dd3;  n1 = nn3;
    
    n23 = n2*n3;
    nm=n1*n23;
    n = nm;
    
    rd1 = 1./d1;
    rd2 = 1./d2;
    rd3 = 1./d3;

    j1 = max(min(xs[0],n1-2),0);
    j2 = max(min(xs[1],n2-2),0);
    j3 = max(min(xs[2],n3-2),0);

    /* should probably check the return value of malloc */
    ttime  = (float *) malloc (n*sizeof(float));
    mask  = (unsigned char *) malloc (n*sizeof(unsigned char));

    slow[0] = 1./slow[1];
    slowv[0] = slowv[1];
    eta[0] = eta[1];
    ttime[0] = SF_HUGE;
    mask[1] = FMM_OUT;

    for (i=1; i<n; i++) {
	ttime[i] = SF_HUGE;
	mask[i] = FMM_OUT;
	slow[i] = 1./slow[i];
	slowv[i] = slowv[i];
    }

    j = AT(j1,j2,j3);
    s2 = slow[j];
    sl = sqrt(s2);
    sv2 = 1./slowv[j];
    svl = sqrt(sv2);
    etal = eta[j];
    sh = sl/sqrt(1+2*etal);
    ttime[j] = 0.;
    mask[j] = FMM_IN;
    nm = n-1;

    /* nh is an estimate of the maximum front size */
    nh = 0;
    if (n1 > 1) nh += 2*n2*n3;
    if (n2 > 1) nh += 2*n1*n3;
    if (n3 > 1) nh += 2*n1*n2;

    sf_pqueue_init(nh);
    sf_pqueue_start();

    /* initialize source */
    for (i1 = j1-1; i1 <= j1+1; i1+=2) {
	if (i1 >= 0 && i1 < n1) {
	    j = AT(i1,j2,j3); INSERT;
	    *pt = sh * d1;
	}
    }
    for (i2 = j2-1; i2 <= j2+1; i2+=2) {
	if (i2 >= 0 && i2 < n2) {
	    j = AT(j1,i2,j3); INSERT;
	    *pt = sh * d2;
	}
    }
    for (i3 = j3-1; i3 <= j3+1; i3+=2) {
	if (i3 >= 0 && i3 < n3) {
	    j = AT(j1,j2,i3); INSERT;
	    *pt = svl * d3; 
	}
    }
    for (i2 = j2-1; i2 <= j2+1; i2+=2) {
	if (i2 >= 0 && i2 < n2) {
	    for (i1 = j1-1; i1 <= j1+1; i1+=2) {
		if (i1 >= 0 && i1 < n1) {
		    j = AT(i1,i2,j3); INSERT;
		    *pt = sh * sqrt(d1*d1+d2*d2);
		}
	    }
	}
    }

    dist2 = d1*d1+d3*d3;
    cost2 = d3*d3/dist2;
    sint2 = 1-cost2;
    ss=2*s2*sv2/(cost2*s2+sint2*(sv2+2*etal*sv2) +
		 sqrt(-8*cost2*etal*s2*sint2*sv2+(cost2*s2+(1 +2*etal)*sint2*sv2)
		      *(cost2*s2+(1+2*etal)*sint2*sv2)));
    dist2 = sqrt(ss*dist2);
    for (i3 = j3-1; i3 <= j3+1; i3+=2) {
	if (i3 >= 0 && i3 < n3) {
	    for (i1 = j1-1; i1 <= j1+1; i1+=2) {
		if (i1 >= 0 && i1 < n1) {
		    j = AT(i1,j2,i3); INSERT;
		    *pt = dist2;
		}
	    }
	}
    }

    dist2 = d2*d2+d3*d3;
    cost2 = d3*d3/dist2;
    sint2 = 1-cost2;
    ss=2*s2*sv2/(cost2*s2+sint2*(sv2+2*etal*sv2) +
		 sqrt(-8*cost2*etal*s2*sint2*sv2+(cost2*s2+(1 +2*etal)*sint2*sv2)
		      *(cost2*s2+(1+2*etal)*sint2*sv2)));
    dist2 = sqrt(ss*dist2);
    for (i3 = j3-1; i3 <= j3+1; i3+=2) {
	if (i3 >= 0 && i3 < n3) {
	    for (i2 = j2-1; i2 <= j2+1; i2+=2) {
		if (i2 >= 0 && i2 < n2) {
		    j = AT(j1,i2,i3); INSERT;
		    *pt = dist2;
		}
	    }
	}
    }

    dist2 = d1*d1+d2*d2+d3*d3;
    cost2 = d3*d3/dist2;
    sint2 = 1-cost2;
    ss=2*s2*sv2/(cost2*s2+sint2*(sv2+2*etal*sv2) +
		 sqrt(-8*cost2*etal*s2*sint2*sv2+(cost2*s2+(1 +2*etal)*sint2*sv2)
		      *(cost2*s2+(1+2*etal)*sint2*sv2)));
    dist2 = sqrt(ss*dist2);
    for (i3 = j3-1; i3 <= j3+1; i3+=2) {
	if (i3 >= 0 && i3 < n3) {
	    for (i2 = j2-1; i2 <= j2+1; i2+=2) {
		if (i2 >= 0 && i2 < n2) {
		    for (i1 = j1-1; i1 <= j1+1; i1+=2) {
			if (i1 >= 0 && i1 < n1) {
			    j = AT(i1,i2,i3); INSERT;
			    *pt = dist2;
			}
		    }
		}
	    }
	}
    }
    /* source initialized */
    sf_warning("v=%f ss=%f dist2=%f nm=%d",slow[1],ss,dist2,nm);

    /* precompute some of the finite-difference coefficients */
    dd[0] = 0.;
    dd[1] = d1*d1;
    dd[2] = d2*d2;
  
    d1 = 1./dd[1];
    d2 = 1./dd[2];
    d3 = 1./(d3*d3);

    dd[3] = 1./(d1+d2);

    if (order==2) 
    {
	while (nm > 0) {
	    pt = sf_pqueue_extract ();
	    i = pt-ttime;
	    i3 = i%n3;
	    i2 = (i/n3)%n2;
	    i1 = i/n23;
	    sf_warning("i1=%d i2=%d i3=%d i=%d time=%f slow=%f slow=%f eta=%f nm=%d",i1,i2,i3,i,ttime[i],slow[i],slowv[i],eta[i],nm);
      
	    *(pm = mask+i) = FMM_IN;
	    if (i1 < n1-1 && *(pm+n23) != FMM_IN) update2 (i1+1,i2,i3, pt+n23, pm+n23, slow[i+n23],slowv[i+n23],eta[i+n23]); 
	    if (i1 > 0    && *(pm-n23) != FMM_IN) update2 (i1-1,i2,i3, pt-n23, pm-n23, slow[i-n23],slowv[i-n23],eta[i-n23]);
	    if (i2 < n2-1 && *(pm+ n3) != FMM_IN) update2 (i1,i2+1,i3, pt+ n3, pm+ n3, slow[i+ n3],slowv[i+ n3],eta[i+ n3]); 
	    if (i2 > 0    && *(pm- n3) != FMM_IN) update2 (i1,i2-1,i3, pt- n3, pm- n3, slow[i- n3],slowv[i- n3],eta[i- n3]);
	    if (i3 < n3-1 && *(pm+  1) != FMM_IN) update2 (i1,i2,i3+1, pt+  1, pm+  1, slow[i+  1],slowv[i+  1],eta[i+  1]); 
	    if (i3 > 0    && *(pm-  1) != FMM_IN) update2 (i1,i2,i3-1, pt-  1, pm-  1, slow[i-  1],slowv[i-  1],eta[i-  1]);
	}
    }

    for (i=0; i<n; i++) {
        time[i] = ttime[i];
	}
    
    
    sf_pqueue_close ();
    free (mask);
}


///////////////////////////////////////////////////////////////////////////


void updatelll (int p1, int p2, int p3, float* tj, char* mj, float s, float rsv, float eta)
{
    float r, b, c, t1=0., t2=0., u, d3r= d3*s*rsv, b2;
    float d3rr,c1,b1,f,g,gg,ggg,gggg,ggggg,ff,f1;
    float A0,A1,A2,m1,m2,m3;
    double tp1,tp2,den;
    unsigned int k, i;

    dd[4] = 1./d3r;

    dd[5] = 1./(d1+d3r);
    dd[6] = 1./(d2+d3r);

    dd[7] = 1./(d1+d2+d3r);

    b = c = 0.; i = k = 0; f = 0.;
    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.;
    tp1 = 0.; tp2=0.;
    if ((p1 > 0   ) && *(mj-n23) && ((t1 = *(tj-n23)) < *tj)) i |= 0x01;
    if ((p1 < n1-1) && *(mj+n23) && ((t2 = *(tj+n23)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;} 
    if (i & 0x01) {
	u = d1*t1; b += u; c += u*t1; f += d1;
	i ^= 0x01; k |= 0x01;
    }
    if ((p2 > 0   ) && *(mj-n3) && ((t1 = *(tj-n3)) < *tj)) i |= 0x01;
    if ((p2 < n2-1) && *(mj+n3) && ((t2 = *(tj+n3)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;} 
    if (i & 0x01) {
	u = d2*t1; b += u; c += u*t1; f += d2;
	i ^= 0x01; k |= 0x02;
    }

    if(k){
	a = dd[k];
	r = b*a;
	tp = r + sqrt(SF_MAX(s*a + r*r - c*a,0));
	den = tp*f-b;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	tp1 = -(c-tp*(2.*b-tp*f))*den;
	tp2 = (tp1*(2.*(b-tp*f)-0.5*tp1*f))*den;
    }

    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.; b2 = 0.; f1 = 0.;
    d3rr = -d3*c*rsv;
    if ((p3 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) i |= 0x01;
    if ((p3 < n3-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;}  
    if (i & 0x01) {
	if ((i & 0x02) && (!(i & 0x01) || t2 > t1)) t1 = t2;
	u = d3rr*t1; b1 = b+u; c1= c+u*t1; g= b+t1*f; gg = f; ggg=t1*t1*b;
	gggg = t1*t1*f; ggggg = t1*b; ff =f+d3r;
	u = d3r*t1; b += u; c += u*t1; f += d3rr; b2 = b; f1 = f;
	i ^= 0x01; k |= 0x04;

	a = dd[k];
	r = b*a;
	tp = r + sqrt(SF_MAX(s*a + r*r - c*a,0));
	den = tp*ff-b2;
	/*sf_warning("tp=%f aaa=%f bbb=%f d3r=%f d1=%f d2=%f d3=%f",tp,ff,b,d3r,d1,d2,d3);*/
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	m1 = rsv*d3*(gggg+4.*ggggg);
	m2 = f1-m1+tp*d3*rsv*(3.*g-2.*tp*gg);
	m3 = b1-tp*m2-rsv*d3*ggg;
	tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(-(c1-2.*tp*b1+tp*tp*f1+tp*tp*tp*d3*rsv*
					  (2*g-tp*gg)+rsv*d3*tp*(2*ggg-tp*gggg-4.*tp*ggggg))*den),SF_ABS(tp1));
	tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS((2.*tp1*(b1-tp*f1-rsv*d3*(ggg-tp*gggg-4.*tp*ggggg+
								  tp*tp*(3.*g-2.*tp*gg)))-0.5*tp1*tp1*ff)*den),SF_ABS(tp2));
    }
    if (!k) return;
    /*tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(tp1),.1);
      tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS(tp2),1.);
      tp3 = SF_SIG(tp3)*SF_MIN(SF_ABS(tp3),1.);*/
    A0=tp; A1=A0+tp1*eta; A2=A1+tp2*eta*eta;
    /*sf_warning("tp=%f tp1=%f tp2=%f A1=%f A2=%f A3=%f den=%f",tp,tp1,tp2,A1,A2,A3,den);*/
    if(SF_ABS(tp1)>0.0000001 && SF_ABS(tp2)>0.0000001){
	den = tp1-tp2*eta;
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	tp = (tp*tp1-(tp*tp2-tp1*tp1)*eta)*den;
    }

    /*sf_warning("tp=%f den=%f k=%d j1=%d j2=%d j3=%d",tp,den,k,p1,p2,p3);*/
    /*  tp = r + sqrt (fabs(tp*tp + r*r - c*a)); */
    if (tp < *tj) {
	*tj = tp;
	if (*mj == FMM_OUT) {
	    nm--; 
	    *mj = FMM_FRONT; 
	    sf_pqueue_insert (tj);
	} 
    }
}

void updateold (int p1, int p2, int p3, float* tj, char* mj, float s, float rsv, float eta)
{
    float r, b, c, t1=0., t2=0., u, d3r= d3*s*rsv, b2;
    float d3rr,c1,b1,f,g,gg,ggg,gggg,ggggg,ff,f1;
    float m1,m2,m3;
    double tp1,tp2,den;
    unsigned int k, i;

    dd[4] = 1./d3r;

    dd[5] = 1./(d1+d3r);
    dd[6] = 1./(d2+d3r);

    dd[7] = 1./(d1+d2+d3r);

    b = c = 0.; i = k = 0; f = 0.;
    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.;
    tp1 = 0.; tp2=0.;
    if ((p1 > 0   ) && *(mj-n23) && ((t1 = *(tj-n23)) < *tj)) i |= 0x01;
    if ((p1 < n1-1) && *(mj+n23) && ((t2 = *(tj+n23)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;} 
    if (i & 0x01) {
	u = d1*t1; b += u; c += u*t1; f += d1;
	i ^= 0x01; k |= 0x01;
    }
    if ((p2 > 0   ) && *(mj-n3) && ((t1 = *(tj-n3)) < *tj)) i |= 0x01;
    if ((p2 < n2-1) && *(mj+n3) && ((t2 = *(tj+n3)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;} 
    if (i & 0x01) {
	u = d2*t1; b += u; c += u*t1; f += d2;
	i ^= 0x01; k |= 0x02;
    }

    if(k){
	a = dd[k];
	r = b*a;
	tp = r + sqrt(SF_MAX(s*a + r*r - c*a,0));
	den = tp*f-b;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	tp1 = -(c-tp*(2.*b-tp*f))*den;
	tp2 = (tp1*(2.*(b-tp*f)-0.5*tp1*f))*den;
    }

    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.; b2 = 0.; f1 = 0.;
    d3rr = -d3*c*rsv;
    if ((p3 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) i |= 0x01;
    if ((p3 < n3-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;}  
    if (i & 0x01) {
	if ((i & 0x02) && (!(i & 0x01) || t2 > t1)) t1 = t2;
	u = d3rr*t1; b1 = b+u; c1= c+u*t1; g= b+t1*f; gg = f; ggg=t1*t1*b;
	gggg = t1*t1*f; ggggg = t1*b; ff =f+d3r;
	u = d3r*t1; b += u; c += u*t1; f += d3rr; b2 = b; f1 = f;
	i ^= 0x01; k |= 0x04;

	a = dd[k];
	r = b*a;
	tp = r + sqrt(SF_MAX(s*a + r*r - c*a,0));
	den = tp*ff-b2;
	/*sf_warning("tp=%f aaa=%f bbb=%f d3r=%f d1=%f d2=%f d3=%f",tp,ff,b,d3r,d1,d2,d3);*/
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	m1 = rsv*d3*(gggg+4.*ggggg);
	m2 = f1-m1+tp*d3*rsv*(3.*g-2.*tp*gg);
	m3 = b1-tp*m2-rsv*d3*ggg;
	tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(-(c1-2.*tp*b1+tp*tp*f1+tp*tp*tp*d3*rsv*
					  (2*g-tp*gg)+rsv*d3*tp*(2*ggg-tp*gggg-4.*tp*ggggg))*den),SF_ABS(tp1));
	tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS((2.*tp1*(b1-tp*f1-rsv*d3*(ggg-tp*gggg-4.*tp*ggggg+
								  tp*tp*(3.*g-2.*tp*gg)))-0.5*tp1*tp1*ff)*den),SF_ABS(tp2));
    }
    if (!k) return;

    if(SF_ABS(tp1)>0.0000001 && SF_ABS(tp2)>0.0000001){
	den = tp1-tp2*eta;
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	tp += tp1*tp1*eta*den;
    }

    if (tp < *tj) {
	*tj = tp;
	if (*mj == FMM_OUT) {
	    nm--; 
	    *mj = FMM_FRONT; 
	    sf_pqueue_insert (tj);
	} 
    }
}

void updateh (int p1, int p2, int p3, float* tj, char* mj, float s, float rsv, float eta)
{
    float r, b, c, t1=0., t2=0., u, d3r= d3*s*rsv, b2;
    float d3rr,c1,b1,f,g,gg,ggg,gggg,ggggg,ff,f1;
    float m1,m2,m3,tpA;
    double tp1,tp2,tp3=0.,tp4=0.,den;
    unsigned int k, i;

    dd[4] = 1./d3r;

    dd[5] = 1./(d1+d3r);
    dd[6] = 1./(d2+d3r);

    dd[7] = 1./(d1+d2+d3r);

    b = c = 0.; i = k = 0; f = 0.;
    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.;
    tp1 = 0.; tp2=0.;
    if ((p1 > 0   ) && *(mj-n23) && ((t1 = *(tj-n23)) < *tj)) i |= 0x01;
    if ((p1 < n1-1) && *(mj+n23) && ((t2 = *(tj+n23)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;} 
    if (i & 0x01) {
	u = d1*t1; b += u; c += u*t1; f += d1;
	i ^= 0x01; k |= 0x01;
    }
    if ((p2 > 0   ) && *(mj-n3) && ((t1 = *(tj-n3)) < *tj)) i |= 0x01;
    if ((p2 < n2-1) && *(mj+n3) && ((t2 = *(tj+n3)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;} 
    if (i & 0x01) {
	u = d2*t1; b += u; c += u*t1; f += d2;
	i ^= 0x01; k |= 0x02;
    }

    if(k){
	a = dd[k];
	r = b*a;
	tp = r + sqrt(SF_MAX(s*a + r*r - c*a,0));
	den = tp*f-b;
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	tp1 = -(c-tp*(2.*b-tp*f))*den;
	tp2 = (tp1*(2.*(b-tp*f)-0.5*tp1*f))*den;
	tp3 = -(tp1*tp1*f-tp2*(2.*(b-tp*f)-tp1*f))*den;
	tp4 = -(tp2*f*(2*tp1+0.5*tp2)-tp3*(2.*(b-tp*f)-tp1*f))*den;
    }

    g = gg = ggg = gggg = ggggg = ff = b1 = c1 = 0.; b2 = 0.; f1 = 0.;
    d3rr = -d3*c*rsv;
    if ((p3 > 0   ) && *(mj-1) && ((t1 = *(tj-1)) < *tj)) i |= 0x01;
    if ((p3 < n3-1) && *(mj+1) && ((t2 = *(tj+1)) < *tj) && 
	((i ^ 0x01) || t2 > t1)) {i |= 0x01; t1 = t2;}  
    if (i & 0x01) {
	if ((i & 0x02) && (!(i & 0x01) || t2 > t1)) t1 = t2;
	u = d3rr*t1; b1 = b+u; c1= c+u*t1; g= b+t1*f; gg = f; ggg=t1*t1*b;
	gggg = t1*t1*f; ggggg = t1*b; ff =f+d3r;
	u = d3r*t1; b += u; c += u*t1; f += d3rr; b2 = b; f1 = f;
	i ^= 0x01; k |= 0x04;

	a = dd[k];
	r = b*a;
	tp = r + sqrt(SF_MAX(s*a + r*r - c*a,0));
	den = tp*ff-b2;
	/*sf_warning("tp=%f aaa=%f bbb=%f d3r=%f d1=%f d2=%f d3=%f",tp,ff,b,d3r,d1,d2,d3);*/
	den = (SF_ABS(den)<0.0000001 ? SF_SIG(den)*10000000. : 1./den);
	m1 = rsv*d3*(gggg+4.*ggggg);
	m2 = f1-m1+tp*d3*rsv*(3.*g-2.*tp*gg);
	m3 = b1-tp*m2-rsv*d3*ggg;
	tp1 = SF_SIG(tp1)*SF_MIN(SF_ABS(-(c1-2.*tp*b1+tp*tp*f1+tp*tp*tp*d3*rsv*
					  (2*g-tp*gg)+rsv*d3*tp*(2*ggg-tp*gggg-4.*tp*ggggg))*den),SF_ABS(tp1));
	tp2 = SF_SIG(tp2)*SF_MIN(SF_ABS((2.*tp1*(b1-tp*f1-rsv*d3*(ggg-tp*gggg-4.*tp*ggggg+
								  tp*tp*(3.*g-2.*tp*gg)))-0.5*tp1*tp1*ff)*den),SF_ABS(tp2));
	tp3 = SF_SIG(tp3)*SF_MIN(SF_ABS(-(tp1*tp1*(f1-m1+6*rsv*d3*tp*
						   (g-tp*gg))-2.*tp2*m3+tp1*tp2*ff)*den),SF_ABS(tp3));
	tp4 = SF_SIG(tp4)*SF_MIN(SF_ABS(-(2.*tp1*(tp1*tp1*rsv*d3*(g-2.*tp*gg)+tp2*(f1-m1+6.*rsv*tp*d3*(g-tp*gg)))+
					  tp2*tp2*0.5*ff-2.*tp3*m3+tp1*tp3*ff)*den),SF_ABS(tp4));
    }
    if (!k) return;

    if(SF_ABS(tp1)>0.0000001 && SF_ABS(tp2)>0.0000001 && SF_ABS(tp3)>0.0000001 && SF_ABS(tp4)>0.0000001){
	den = tp1-tp2*eta;
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	tpA = tp+tp1*tp1*eta*den;
	den = (tp2-eta*tp3)*(eta*tp1*tp3*tp3*tp3+eta*eta*tp2*tp2*tp3*tp4+ 
			     tp2*tp2*tp2*(-tp3+eta*tp4)+tp2*tp3*(tp1*tp3-eta*eta*tp3*tp3-2*eta*tp1*tp4));
	den = (SF_ABS(den)<0.000000001 ? SF_SIG(den)*1000000000. : 1./den);
	tp += eta*(eta*tp2*tp2*tp2*tp2*tp2*(-tp3 +eta*tp4)+ 
		   tp1*tp2*tp2*tp2*(tp2-2*eta*tp3)*(-tp3+eta*tp4)+ 
		   tp1*tp1*tp3*(tp2-eta*tp3)*(tp2*tp3+eta*tp3*tp3-2*eta*tp2*tp4))*den;
	sf_warning("tp=%f tpA=%f",tp,tpA);
    }

    if (tp < *tj) {
	*tj = tp;
	if (*mj == FMM_OUT) {
	    nm--; 
	    *mj = FMM_FRONT; 
	    sf_pqueue_insert (tj);
	} 
    }
}















 




///////////////////////////////////////////////////////////////////////////
///////////////      INITIAL-VALUE RAY TRACER     /////////////////////////
///////////////////////////////////////////////////////////////////////////


void FSM3DTTI(GridClass G, double *shot, double *Vp, double *Eta, double *out)
{
	
    int n1, n2, n3, b1, b2, b3; 
    int i, j, k; 
    int nshot, ndim;
    int is,order,n123, *p;
    int order=2;                    // [1,2] Accuracy order

    
// int sy, sz, loc=0; y
    //int i, j, loop, ndim, niter, nfpi, fac;
    float o1, o2, o3, d1, d2, d3;
    float br1, br2, br3;
    float *s, *t, *eta, *v, *vv, *q;
//    float *s, *t, *tau, *vz, *vx, *eps, *del, *theta, *T0, *py0, *pz0, *tn, *tn1, a0, b0, c0, vx0, vz0, st0, ct0;
//    float *py, *pz, pydash, pzdash, *rhs, sum;
    //bool optloc;
    bool plane[3];
   
    n1 = G.nz;  n2 = G.nx;  n3 = G.ny;
 
    d1 = G.dz;  d2 = G.dx;  d3 = G.dy;
 
	o1 = G.z0;  o2 = G.x0;  o3 = G.y0;

    //niter = 4;       /* number of sweeping iterations */
 
    //nfpi = 3;        /* number of fixed-point iterations */

    //fac = 1;         /* Type of factorization: (0)Additive, (1)Multiplicative */
                     /* Multiplicative factorization is more stable */
    //optloc = false;  /* Selects optimal location for homogeneous medium parameter */
                     /* Useful for stability of additive factorization when the highest velocity in the medium is 
                        much larger than the velocity at the source point */
    
    n123 = n1*n2*n3;
    
//    vvf = 
//    eta = 
            


    /* Constant-velocity box around the source (in physical dimensions) */
    br1=d1;     br2=d2;     br3=d3;

    /* plane-wave source */
    plane[2]=false;
    plane[1]=false;
    plane[0]=false;
    

    /* Constant-velocity box around the source (in samples) */
    b1= plane[2]? n1: (int) (br1/d1+0.5); 
    b2= plane[1]? n2: (int) (br2/d2+0.5); 
    b3= plane[0]? n3: (int) (br3/d3+0.5); 


    if( b1<1 ) b1=1;  
    if( b2<1 ) b2=1;  
    if( b3<1 ) b3=1;

    nshot = 1;
	ndim = 3;
    

    

    s[0][0]=0.; //z
    s[0][1]=o2 + 0.5*(n2-1)*d2; //y
    s[0][2]=o3 + 0.5*(n3-1)*d3; //x

    t  = new float[n123];
    v  = new float[n123];           //NMP velocity squared
    vv = new float[n123];           // vertical velocity squared
    q  = new float[n123];           // eta
    p  = new int[n123];
    
    /* Reading input parameters */
    /* input has size G.nx x G.ny x G.nz */
    
    for (i=0;i<G.nx;i++){
        for (j=0;j<G.ny;j++){
            for (k=0;k<G.nz;k++){
                v[k+j*G.nz+i*G.nz*G.ny] = (float)Vp[i+j*G.nx+k*G.nx*G.ny];
                q[k+j*G.nz+i*G.nz*G.ny] = (float)Eta[i+j*G.nx+k*G.nx*G.ny];
                vv
            }
        }
    }
    
//    for(i = 0; i < n123; i++) {
//    	    v[i] = v[i]*v[i];
//	  vv[i] = vv[i]*vv[i];


      fastvti(t,v,vv,q,p, plane,
		  n3,n2,n1,
		  o3,o2,o1,
		  d3,d2,d1,
		  shot[0],shot[1],shot[2], 
		  b3,b2,b1,
		  order);
       
    
    /* Convert angles from degrees to radians */
    //    vx[i] = vz[i]*sqrtf(1+2*eps[i]); /* Compute horizontal velocity */
    /* Converting source location into grid points */
    //sy = (int)((shot[0] - o2) / d2 + 0.5f);
	//sz = (int)((shot[2] - o1) / d1 + 0.5f);
    //mexErrMsgTxt("Choose fac=0(Additive) or fac=1(Multiplicative) factorization");
    
    
    //Fill output array
    for (i=0;i<G.nx;i++){
        for (j=0;j<G.ny;j++){
            for (k=0;k<G.nz;k++){
                out[i+j*G.nx+k*G.nx*G.ny] = (double)t[k+j*G.nz+i*G.nz*G.ny];
            }
        }
    }
    
    delete[] t, v, vv, q, p; 

    return;
}
    
// ///////////////////////////////////////////////////////////////////////////
// //////////////////      Gateway routine      //////////////////////////////
// ///////////////////////////////////////////////////////////////////////////

// Transfer Data and results from MATLAB to C++ and back

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

    //declare variables
    const mwSize *dims;
    mwSize ndim;
    double *tmp;
    double *shot, *Vp, *Epsilon, *Delta, *Theta, *tti;
    
    GridClass G;
    
    
    //check input
    if (nrhs!=6) mexErrMsgTxt("Three input arguments are required.");
    if (nlhs!=1) mexErrMsgTxt("One output argument is required.");

     
    // Get G - parameters of the velocity grid
    if (!mxIsClass(prhs[0], "GridClass")) mexErrMsgTxt("Input (1st arg.) must be an object of class GridClass");
    if (mxGetProperty(prhs[0],0,"nx")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'nx' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"nx"));
	G.nx = int(*tmp);
    if (mxGetProperty(prhs[0],0,"ny")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'ny' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"ny"));
	G.ny = int(*tmp);
    if (mxGetProperty(prhs[0],0,"nz")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'nz' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"nz"));
	G.nz = int(*tmp);
    if (mxGetProperty(prhs[0],0,"nt")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'nt' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"nt"));
	G.nt = int(*tmp);
    if (mxGetProperty(prhs[0],0,"x0")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'x0' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"x0"));
	G.x0 = double(*tmp);
    if (mxGetProperty(prhs[0],0,"y0")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'y0' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"y0"));
	G.y0 = double(*tmp);
    if (mxGetProperty(prhs[0],0,"z0")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'z0' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"z0"));
	G.z0 = double(*tmp);
    if (mxGetProperty(prhs[0],0,"t0")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 't0' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"t0"));
	G.t0 = double(*tmp);
    if (mxGetProperty(prhs[0],0,"dx")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'dx' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"dx"));
	G.dx = double(*tmp);
    if (mxGetProperty(prhs[0],0,"dy")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'dy' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"dy"));
	G.dy = double(*tmp);
    if (mxGetProperty(prhs[0],0,"dz")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'dz' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"dz"));
	G.dz = double(*tmp);
    if (mxGetProperty(prhs[0],0,"dt")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'dt' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"dt"));
	G.dt = double(*tmp);
    if (mxGetProperty(prhs[0],0,"mx")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'mx' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"mx"));
	G.mx = double(*tmp);
    if (mxGetProperty(prhs[0],0,"my")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'my' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"my"));
	G.my = double(*tmp);
    if (mxGetProperty(prhs[0],0,"mz")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'mz' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"mz"));
	G.mz = double(*tmp);
    if (mxGetProperty(prhs[0],0,"mt")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'mt' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"mt"));
	G.mt = double(*tmp);

    
    //associate input pointers
    shot = (double *)mxGetPr(prhs[1]);
    Vp = (double *)mxGetPr(prhs[2]);
    Epsilon = (double *)mxGetPr(prhs[3]);
    Delta = (double *)mxGetPr(prhs[4]);
    Theta = (double *)mxGetPr(prhs[5]);
    
    //check if dimensions velmod, Epsilon, Delta and Theta are correct
    ndim = mxGetNumberOfDimensions(prhs[2]);
    dims = mxGetDimensions(prhs[2]);
    if (ndim!=2) mexErrMsgTxt("Number of dimensions of Vp array is not correct.");
    if (dims[0]!=G.nx) mexErrMsgTxt("Vp should have size [G.nx, G.nz].");
    if (dims[1]!=G.nz) mexErrMsgTxt("Vp should have size [G.nx, G.nz].");

    ndim = mxGetNumberOfDimensions(prhs[3]);
    dims = mxGetDimensions(prhs[3]);
    if (ndim!=2) mexErrMsgTxt("Number of dimensions of Epsilon array is not correct.");
    if (dims[0]!=G.nx) mexErrMsgTxt("Epsilon should have size [G.nx, G.nz].");
    if (dims[1]!=G.nz) mexErrMsgTxt("Epsilon should have size [G.nx, G.nz].");

    ndim = mxGetNumberOfDimensions(prhs[4]);
    dims = mxGetDimensions(prhs[4]);
    if (ndim!=2) mexErrMsgTxt("Number of dimensions of Delta array is not correct.");
    if (dims[0]!=G.nx) mexErrMsgTxt("Delta should have size [G.nx, G.nz].");
    if (dims[1]!=G.nz) mexErrMsgTxt("Delta should have size [G.nx, G.nz].");
    
    ndim = mxGetNumberOfDimensions(prhs[5]);
    dims = mxGetDimensions(prhs[5]);
    if (ndim!=2) mexErrMsgTxt("Number of dimensions of Theta array is not correct.");
    if (dims[0]!=G.nx) mexErrMsgTxt("Theta should have size [G.nx, G.nz].");
    if (dims[1]!=G.nz) mexErrMsgTxt("Theta should have size [G.nx, G.nz].");
    
    //create array for output
    plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
    
    //associate output pointer
    tti = (double *)mxGetPr(plhs[0]);
    
    // Call the master routine
    if(info!= 0) mexPrintf("Calling FSM routine from Mex gateway...\n");
    FSM2DTTI(G, shot, Vp, Epsilon, Delta, Theta, tti);
    if(info!= 0) mexPrintf("FSM finished.\n");
}










