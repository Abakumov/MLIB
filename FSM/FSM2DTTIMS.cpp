// Fast sweeping factored TTI eikonal solver (2D)
// Original C++ code by Umair bin Waheed
// Mex wrapper by Ivan Abakumov
// Freie Universität Berlin
// e-mail: abakumov_ivan@mail.ru
// 20th of March 2018
// version: 0.1
//
//////////////////////////////////////////////////////
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

void sf_derivative_2D (float *t,float *dtdz, float *dtdx,
                       double d2, double d1,
                       int n2, int n1, int accuracy)

{
    int i1, i2;

    // For second order accuracy
    if (accuracy==2){

        if (n1<3 || n2<3){
            mexErrMsgTxt("Need more points for second order accuracy");
        }

        for(i2 = 0; i2 < n2; i2++){
            for(i1 = 1; i1 < n1-1; i1++){

            dtdz[n1*i2 + i1] = (t[n1*i2+(i1+1)]-t[n1*i2+(i1-1)])/(2*d1);

            }

        dtdz[n1*i2+0] = (-t[n1*i2 + 2] +4*t[n1*i2 + 1] -3*t[n1*i2 + 0])/(2*d1);

        dtdz[n1*i2+(n1-1)] = (3*t[n1*i2 + n1-1]- 4*t[n1*i2 + n1-2] + t[n1*i2 + n1-3])/(2*d1);

        }


        for(i1 = 0; i1 < n1; i1++){
            for(i2 = 1; i2 < n2-1; i2++){

            dtdx[n1*i2 + i1] = (t[n1*(i2+1)+i1] - t[n1*(i2-1)+i1])/(2*d2);

            }

        dtdx[n1*0+i1] = (-t[n1*2 + i1] +4*t[n1*1 + i1] -3*t[n1*0 + i1])/(2*d2);

        dtdx[n1*(n2-1)+i1] = (3*t[n1*(n2-1) + i1]- 4*t[n1*(n2-2) + i1] + t[n1*(n2-3) + i1])/(2*d2);

        }

    }


    //For fourth order accuracy
    if (accuracy==4){


        if (n1<5 || n2<5){
        mexErrMsgTxt("Need more points for fourth order accuracy");
        return;
        }
        for(i2 = 0; i2 < n2; i2++){
            for(i1 = 2; i1 < n1-2; i1++){

            dtdz[n1*i2 + i1] = (-t[n1*i2+(i1+2)] + 8*t[n1*i2+(i1+1)]-8*t[n1*i2+(i1-1)]+t[n1*i2+(i1-2)])/(12*d1);

            }

        dtdz[n1*i2+0] = (-3*t[n1*i2 + 4] + 16* t[n1*i2 + 3] -36* t[n1*i2 + 2] +48*t[n1*i2 + 1] -25*t[n1*i2 + 0])/(12*d1);

        dtdz[n1*i2+1] = (t[n1*i2 + 4] - 6* t[n1*i2 + 3] +18* t[n1*i2 + 2]-10*t[n1*i2 + 1] -3*t[n1*i2 + 0])/(12*d1);

        dtdz[n1*i2+(n1-1)] = (25*t[n1*i2 + n1-1]- 48*t[n1*i2 + n1-2] + 36*t[n1*i2 + n1-3] - 16*t[n1*i2 + n1-4] + 3*t[n1*i2 + n1-5])/(12*d1);

        dtdz[n1*i2+(n1-2)] = (3*t[n1*i2 + n1-1] + 10*t[n1*i2 + n1-2]- 18*t[n1*i2 + n1-3] + 6*t[n1*i2 + n1-4] - t[n1*i2 + n1-5])/(12*d1);

        }


        for(i1 = 0; i1 < n1; i1++){
            for(i2 = 2; i2 < n2-2; i2++){

            dtdx[n1*i2 + i1] = (-t[n1*(i2+2)+i1] + 8* t[n1*(i2+1)+i1] -8*t[n1*(i2-1)+i1] + t[n1*(i2-2)+i1])/(12*d2);

            }

        dtdx[n1*0+i1] = (-3*t[n1*4 + i1] +16*t[n1*3 + i1] -36*t[n1*2 + i1] +48*t[n1*1 + i1] -25*t[n1*0 + i1])/(12*d2);

        dtdx[n1*1+i1] = (t[n1*4 + i1] - 6*t[n1*3 + i1] + 18*t[n1*2 + i1] -10*t[n1*1 + i1] -3*t[n1*0 + i1])/(12*d2);


        dtdx[n1*(n2-1)+i1] = (25*t[n1*(n2-1) + i1]- 48*t[n1*(n2-2) + i1] + 36*t[n1*(n2-3) + i1]-16*t[n1*(n2-4) + i1]+3*t[n1*(n2-5) + i1])/(12*d2);

        dtdx[n1*(n2-2)+i1] = (3*t[n1*(n2-1) + i1] + 10*t[n1*(n2-2) + i1]- 18*t[n1*(n2-3) + i1] + 6*t[n1*(n2-4) + i1]- t[n1*(n2-5) + i1])/(12*d2);


        }

    }

   
}

///////////////////////////////////////////////////////////////////////////

bool sf_init_fast_sweep (float *tau,
                         int n2, int n1,
                         float o2, float o1,
                         float d2, float d1,
                         int shoty, int shotz, int fac)
/*< initialize >*/
{
    int i, n12;

    if (NULL == tau)
        return false;

    if (shoty < 0 || shoty >= n2 ||
        shotz < 0 || shotz >= n1)
        return false;

    n12 = n1 * n2;

    for (i = 0; i < n12; i++)
        tau[i] = SF_HUGE;

    if(fac==0) tau[shoty * n1 + shotz] = 0.0;
    else if(fac==1) tau[shoty * n1 + shotz] = 1.0;
    else mexErrMsgTxt("Choose fac=0(Additive) or fac=1(Multiplicative) factorization");

    return true;
}


///////////////////////////////////////////////////////////////////////////



bool isCausalRoot(float root, float uz, float uy, 
	              int sz, int sy,
		          float a, float b, float c,
                  float dz, float dy){

	float dtdz, dtdy, vgz, vgy;
  
	dtdz = (root-uz)/(sz*dz);
	dtdy = (root-uy)/(sy*dy);
	vgz = b*dtdz + c*dtdy;
	vgy = a*dtdy + c*dtdz;


	if(vgz*fabs(dtdz)*sz>0.0 && vgy*fabs(dtdy)*sy>0.0) return true;
	else return false;
}


///////////////////////////////////////////////////////////////////////////

static void sf_fast_sweep_2d_stencil_add (float *tau, float *T0,
                                      float *py0, float *pz0,
                                      float *vz, float *vx,
                                      float *theta, float *rhs,
                                      int i, int j, int n1, int n2,
                                      float d1, float d2,
                                      int shotz, int shoty) {

    /* A variable with c at the end means its value of current grid point.
       For example: vzc is vz at (i,j) or the current grid point */

    /* Sign variable is chosen assuming that the information is spread along 
       the positive y and z dimensions */

    float vzc = vz[j * n1 + i], vxc = vx[j * n1 + i], thetac = theta[j*n1+i], rhsc = rhs[j*n1+i];
    float T0c = T0[j * n1 + i];
    float tauij = tau[j * n1 + i], root;
    int sy=0, sz=0;
    float dz = d1, dz2 = dz*dz, dy = d2, dy2 = dy*dy;
    float vzc2 = vzc*vzc, vxc2 = vxc*vxc;
    float ct = cos(thetac), st = sin(thetac);
    float ct2 = ct*ct, st2 = st*st;    
    float a = vxc2*ct2 + vzc2*st2, b = vxc2*st2 + vzc2*ct2, c = (vxc2-vzc2)*st*ct; 
    float py0c = py0[j * n1 + i], pz0c=pz0[j * n1 + i];
    float tauy, tauz, T0y, T0z, tauC, tauCy, tauCz;
    float num, den, t1, t2;

    if(i==shotz && j==shoty) return;

    if(i==0 || i==n1-1){
        if(i==0){
            tauz = tau[j * n1 + 1]; T0z = T0[j * n1 + 1]; sz = -1; 
        }
        else{
            tauz = tau[j * n1 + n1-2]; T0z = T0[j * n1 + n1-2]; sz = 1; 
        }
    }
    else{
        if((tau[j*n1+i-1]+T0[j*n1+i-1]) < (tau[j*n1+i+1]+T0[j*n1+i+1])){
            tauz = tau[j * n1 + i - 1]; T0z = T0[j * n1 + i - 1]; sz = 1;
        }
		else{
            tauz = tau[j * n1 + i + 1]; T0z = T0[j * n1 + i + 1]; sz = -1; 
        }

    }

    if(j==0 || j==n2-1){
        if(j==0){
            tauy = tau[n1 + i]; T0y = T0[n1 + i]; sy = -1; 
        }
        else{
            tauy = tau[(n2-2)*n1 + i]; T0y = T0[(n2-2)*n1 + i]; sy = 1; 
        }
    }
    else{
        if((tau[(j-1)*n1+i]+T0[(j-1)*n1+i]) < (tau[(j+1)*n1+i]+T0[(j+1)*n1+i])){
            tauy = tau[(j-1)*n1+i]; T0y = T0[(j-1)*n1+i]; sy = 1;
        }
		else{
            tauy = tau[(j+1) * n1 + i ]; T0y = T0[(j+1) * n1 + i ]; sy = -1; 
        }

    }
    

    if (SF_HUGE == tauy && SF_HUGE == tauz)
        return;


    /* Took absolute of py0c and pz0c because point C is in the upwind direction */

	tauCy = tauy - dy*fabs(py0c) + dy*sqrtf((b*rhsc)/(a*b-c*c));
	tauCz = tauz - dz*fabs(pz0c) + dz*sqrtf((a*rhsc)/(a*b-c*c));



    /* The condition used below is different than the one in the isotropic eikonal solver 
       and it holds real importance in getting the correct solution. */

    if (tauy == SF_HUGE) {
        tauC = tauCz;
    } else if (tauz == SF_HUGE) {
        tauC = tauCy;
    } else {

        t1 = a*dz2*tauy + dy2*sz*(-dz*(c*py0c + b*pz0c) + b*sz*tauz) 
              + dy*dz*sy*(-a*dz*py0c + c*(-dz*pz0c + sz*(tauy+tauz)));

        t2 = 2*c*dy*dz*rhsc*sy*sz + b*dy2*rhsc + c*c*(dz*pz0c*sy - dy*py0c*sz 
             + sy*sz*tauy -sy*sz*tauz)*(dz*pz0c*sy - dy*py0c*sz + sy*sz*tauy -sy*sz*tauz) 
             + a*(dz2*rhsc - b*(dz*pz0c*sy - dy*py0c*sz + sy*sz*tauy - sy*sz*tauz)
             *(dz*pz0c*sy - dy*py0c*sz + sy*sz*tauy - sy*sz*tauz));
              

        num = t1 + dy2*dz2*sqrt(t2/(dy2*dz2));

        den = a*dz2 + dy*sz*(2*c*dz*sy + b*dy*sz);


        root = num/den;


	    if(isCausalRoot(root+T0c,tauz+T0z,tauy+T0y,sz,sy,a,b,c,dz,dy)) tauC=root;
        else {
            if(tauCy+T0y < tauCz+T0z) tauC=tauCy;
            else tauC = tauCz;
        }

    }

    if (tauC+T0c < tauij+T0c)
        tau[j * n1 + i] = tauC;
}


///////////////////////////////////////////////////////////////////////////

static void sf_fast_sweep_2d_stencil_mul (float *tau, float *T0,
                                      float *py0, float *pz0,
                                      float *vz, float *vx,
                                      float *theta, float *rhs,
                                      int i, int j, int n1, int n2,
                                      float d1, float d2,
                                      int shotz, int shoty) {

    /* A variable with c at the end means its value of current grid point.
       For example: vzc is vz at (i,j) or the current grid point */

    /* Sign variable is chosen assuming that the information is spread along 
       the positive y and z dimensions */

    float vzc = vz[j * n1 + i], vxc = vx[j * n1 + i], thetac = theta[j*n1+i], rhsc = rhs[j*n1+i];
    float T0c = T0[j * n1 + i], T0c2 = T0c*T0c;
    float tauij = tau[j * n1 + i], root;
    int sy=0, sz=0;
    float dz = d1, dz2 = dz*dz, dy = d2, dy2 = dy*dy;
    float vzc2 = vzc*vzc, vxc2 = vxc*vxc;
    float ct = cos(thetac), st = sin(thetac);
    float ct2 = ct*ct, st2 = st*st;    
    float a = vxc2*ct2 + vzc2*st2, b = vxc2*st2 + vzc2*ct2, c = (vxc2-vzc2)*st*ct; 
    float py0c = py0[j * n1 + i], pz0c=pz0[j * n1 + i];
    float py0c2 = py0c*py0c, pz0c2 = pz0c*pz0c;
    float tauy, tauz, T0y, T0z, tauC, tauCy, tauCz;
    float num, den, t1, t2;

    if(i==shotz && j==shoty) return;

    if(i==0 || i==n1-1){
        if(i==0){
            tauz = tau[j * n1 + 1]; T0z = T0[j * n1 + 1]; sz = -1; 
        }
        else{
            tauz = tau[j * n1 + n1-2]; T0z = T0[j * n1 + n1-2]; sz = 1; 
        }
    }
    else{
        if((tau[j*n1+i-1]*T0[j*n1+i-1]) < (tau[j*n1+i+1]*T0[j*n1+i+1])){
            tauz = tau[j * n1 + i - 1]; T0z = T0[j * n1 + i - 1]; sz = 1;
        }
		else{
            tauz = tau[j * n1 + i + 1]; T0z = T0[j * n1 + i + 1]; sz = -1; 
        }

    }

    if(j==0 || j==n2-1){
        if(j==0){
            tauy = tau[n1 + i]; T0y = T0[n1 + i]; sy = -1; 
        }
        else{
            tauy = tau[(n2-2)*n1 + i]; T0y = T0[(n2-2)*n1 + i]; sy = 1; 
        }
    }
    else{
        if((tau[(j-1)*n1+i]*T0[(j-1)*n1+i]) < (tau[(j+1)*n1+i]*T0[(j+1)*n1+i])){
            tauy = tau[(j-1)*n1+i]; T0y = T0[(j-1)*n1+i]; sy = 1;
        }
		else{
            tauy = tau[(j+1) * n1 + i ]; T0y = T0[(j+1) * n1 + i ]; sy = -1; 
        }

    }
    

    if (SF_HUGE == tauy && SF_HUGE == tauz)
        return;


    /* Took absolute of py0c and pz0c because point C is in the upwind direction */

    tauCy = (T0c*tauy + dy*sqrtf(b*rhsc/(a*b-c*c)))/(T0c + fabs(py0c)*dy);
    tauCz = (T0c*tauz + dz*sqrtf(a*rhsc/(a*b-c*c)))/(T0c + fabs(pz0c)*dz);


    // The condition used below is different than the one in the isotropic eikonal solver 
    // and it holds real importance in getting the correct solution.

    if (tauy == SF_HUGE) {
        tauC = tauCz;
    } else if (tauz == SF_HUGE) {
        tauC = tauCy;
    } else {

        t1 = a*dz2*T0c2*tauy + dy2*sz*T0c*(c*dz*py0c + b*dz*pz0c + b*sz*T0c)*tauz 
             + dy*dz*sy*T0c*(a*dz*py0c*tauy + c*(dz*pz0c*tauy + sz*T0c*(tauy+tauz)));

        t2 = 2*c*dy*dz*rhsc*(dy*py0c+sy*T0c)*(dz*pz0c+sz*T0c) 
             + b*dy2*rhsc*(dz*pz0c + sz*T0c)*(dz*pz0c + sz*T0c) + c*c*T0c*T0c*(dz*pz0c*sy*tauy 
             + sy*sz*T0c*(tauy-tauz) - dy*py0c*sz*tauz)*(dz*pz0c*sy*tauy + sy*sz*T0c*(tauy-tauz) 
             - dy*py0c*sz*tauz) + a*(-T0c2*(dz2*(-rhsc+b*pz0c2*tauy*tauy) 
             + 2*b*dz*pz0c*sz*T0c*tauy*(tauy-tauz) + b*T0c2*(tauy-tauz)*(tauy-tauz)) 
             + 2*dy*py0c*sy*T0c*(dz2*rhsc + b*dz*pz0c*sz*T0c*tauy*tauz + b*T0c2*(tauy-tauz)*tauz) 
             + dy2*py0c2*(dz2*rhsc - b*T0c2*tauz*tauz));
           

        num = t1 + dy2*dz2*sqrtf(t2/(dy2*dz2));

        den = a*dz2*(dy*py0c + sy*T0c)*(dy*py0c + sy*T0c) + dy*(dz*pz0c + sz*T0c)*(2*c*dz*(dy*py0c 
              + sy*T0c) + b*dy*(dz*pz0c + sz*T0c));


        root = num/den;



	    if(isCausalRoot(root*T0c,tauz*T0z,tauy*T0y,sz,sy,a,b,c,dz,dy)) tauC=root;
        else {
            if(tauCy*T0y < tauCz*T0z) tauC=tauCy;
            else tauC = tauCz;
        }

    }

    if (tauC*T0c < tauij*T0c)
        tau[j * n1 + i] = tauC;
}


///////////////////////////////////////////////////////////////////////////


void sf_run_fast_sweep (float *tau, float *T0, 
                        float *py0, float *pz0,
                        float *vz, float *vx,
                        float *theta, float *rhs,
                        int niter,
                        int n2, int n1,
                        float o2, float o1,
                        float d2, float d1,
                        int shoty, int shotz, int fac)
/*< run sweeps >*/
{
    int i, j , l = 0;

    if(fac==0){

        for (l = 0; l < niter; l++) {

            for (j = n2 - 1; j >= 0; j--) {
                for (i = 0; i < n1; i++) {
                    sf_fast_sweep_2d_stencil_add (tau, T0, py0, pz0, vz, vx,
                                              theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }
            for (j = n2 - 1; j >= 0; j--) {
                for (i = n1 - 1; i >= 0; i--) {
                    sf_fast_sweep_2d_stencil_add (tau, T0, py0, pz0, vz, vx,
                                              theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }
            for (j = 0; j < n2; j++) {
                for (i = 0; i < n1; i++) {
                    sf_fast_sweep_2d_stencil_add (tau, T0, py0, pz0, vz, vx,
                                          theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }


            for (j = 0; j < n2; j++) {
                for (i = n1 - 1; i >= 0; i--) {
                    sf_fast_sweep_2d_stencil_add (tau, T0, py0, pz0, vz, vx,
                                              theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }

        }

    }
    else if(fac==1){

        for (l = 0; l < niter; l++) {

            for (j = n2 - 1; j >= 0; j--) {
                for (i = 0; i < n1; i++) {
                    sf_fast_sweep_2d_stencil_mul (tau, T0, py0, pz0, vz, vx,
                                              theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }
            for (j = n2 - 1; j >= 0; j--) {
                for (i = n1 - 1; i >= 0; i--) {
                    sf_fast_sweep_2d_stencil_mul (tau, T0, py0, pz0, vz, vx,
                                              theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }
            for (j = 0; j < n2; j++) {
                for (i = 0; i < n1; i++) {
                    sf_fast_sweep_2d_stencil_mul (tau, T0, py0, pz0, vz, vx,
                                          theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }


            for (j = 0; j < n2; j++) {
                for (i = n1 - 1; i >= 0; i--) {
                    sf_fast_sweep_2d_stencil_mul (tau, T0, py0, pz0, vz, vx,
                                              theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }

        }

    }
    else mexErrMsgTxt("Choose fac=0(Additive) or fac=1(Multiplicative) factorization");
}


///////////////////////////////////////////////////////////////////////////
///////////////      INITIAL-VALUE RAY TRACER     /////////////////////////
///////////////////////////////////////////////////////////////////////////


void FSM2DTTI(GridClass G, double *shot, double *Vp, double *Epsilon, double *Delta, double *Theta, double *out)
{
	
    int n1, n2, n3, i, j, loop, n12, ndim, niter, nfpi, fac;
    float o1, o2, o3, d1, d2, d3;
    float *s, *t, *tau, *vz, *vx, *eps, *del, *theta, *T0, *py0, *pz0, *tn, *tn1, a0, b0, c0, vx0, vz0, st0, ct0;
    int sy, sz, loc=0; 
    float *py, *pz, pydash, pzdash, *rhs, sum;
    bool optloc;
   
    n1 = G.nz;  n2 = G.nx;  n3 = G.ny;
 
    d1 = G.dz;  d2 = G.dx;  d3 = G.dy;
 
	o1 = G.z0;  o2 = G.x0;  o3 = G.y0;

    niter = 4;       /* number of sweeping iterations */
 
    nfpi = 3;        /* number of fixed-point iterations */

    fac = 1;         /* Type of factorization: (0)Additive, (1)Multiplicative */
                     /* Multiplicative factorization is more stable */
    optloc = false;  /* Selects optimal location for homogeneous medium parameter */
                     /* Useful for stability of additive factorization when the highest velocity in the medium is 
                        much larger than the velocity at the source point */
    n12 = n1*n2;

    t  = new float[n12];
    tau  = new float[n12];
    vz  = new float[n12];
    vx  = new float[n12];
    eps  = new float[n12];
    del  = new float[n12];
    theta  = new float[n12];
    T0  = new float[n12];
    py0  = new float[n12];
    pz0  = new float[n12];
    py  = new float[n12];
    pz  = new float[n12];
    rhs = new float[n12];
    tn = new float[n12];
    tn1 = new float[n12];
 
    /* Reading input parameters */
    /* input has size G.nx x G.nz */
    
    for (i=0;i<G.nx;i++){
        for (j=0;j<G.nz;j++){
            vz[j+i*G.nz] = (float)Vp[i+j*G.nx];
            eps[j+i*G.nz] = (float)Epsilon[i+j*G.nx];
            del[j+i*G.nz] = (float)Delta[i+j*G.nx];
            theta[j+i*G.nz] = (float)Theta[i+j*G.nx];
        }
    }
    
    /* Convert angles from degrees to radians */
 
    for(i=0;i<n12;i++){
         
        vx[i] = vz[i]*sqrtf(1+2*eps[i]); /* Compute horizontal velocity */
        theta[i] = theta[i]*PI/180.0;       /* Convert angle from degrees to radians */
        rhs[i] = 1.0;     /* Inititialize rhs to 1 */
        tn[i] = 0.;       /* tn is the current time, and tn1 is the time from previous iteration */
 
    }
 
    /* Converting source location into grid points */
    sy = (int)((shot[0] - o2) / d2 + 0.5f);
	sz = (int)((shot[2] - o1) / d1 + 0.5f);
 
    for(loop=0;loop<nfpi;loop++){
     
        if(optloc) loc = (int)(n1/2.);
            
        if(optloc){
            vx0 = vx[sy*n1+loc]; vz0 = vz[sy*n1+loc]; 
            ct0 = cos(theta[sy*n1+loc]); st0 = sin(theta[sy*n1+loc]);
 
            a0 = (vx0*vx0*ct0*ct0 + vz0*vz0*st0*st0)/rhs[sy*n1+loc];
            b0 = (vx0*vx0*st0*st0 + vz0*vz0*ct0*ct0)/rhs[sy*n1+loc];
            c0 = ((vx0*vx0 - vz0*vz0)*st0*ct0)/rhs[sy*n1+loc];
 
        }
        else{
            vx0 = vx[sy*n1+sz]; vz0 = vz[sy*n1+sz]; 
            ct0 = cos(theta[sy*n1+sz]); st0 = sin(theta[sy*n1+sz]);
 
            a0 = (vx0*vx0*ct0*ct0 + vz0*vz0*st0*st0)/rhs[sy*n1+sz];
            b0 = (vx0*vx0*st0*st0 + vz0*vz0*ct0*ct0)/rhs[sy*n1+sz];
            c0 = ((vx0*vx0 - vz0*vz0)*st0*ct0)/rhs[sy*n1+sz];
 
        }

        /* Traveltime and derivative computation for homogeneous TEA medium */
        for(i=0;i<n1;i++){
            for(j=0;j<n2;j++){

                if(i==sz && j==sy){
                    T0[j*n1+i] = 0.; 
                    py0[j*n1+i] =0;  pz0[j*n1+i]=0;
                    continue;
                } 


                T0[j*n1+i] = sqrtf((b0*(j-sy)*(j-sy)*d2*d2 - 2*c0*(j-sy)*(i-sz)*d1*d2 
                           + a0*(i-sz)*(i-sz)*d1*d1)/(a0*b0-c0*c0));
 

                py0[j*n1+i] = (b0*(j-sy)*d2 - c0*(i-sz)*d1)/(sqrtf((b0*(j-sy)*(j-sy)*d2*d2
                            - 2*c0*(j-sy)*(i-sz)*d1*d2 + a0*(i-sz)*(i-sz)*d1*d1)
        	    			   *(a0*b0-c0*c0)));
 
                pz0[j*n1+i] = (a0*(i-sz)*d1 - c0*(j-sy)*d2)/(sqrtf((b0*(j-sy)*(j-sy)*d2*d2
                            - 2*c0*(j-sy)*(i-sz)*d1*d2 + a0*(i-sz)*(i-sz)*d1*d1)
    				    	   *(a0*b0-c0*c0)));

            }
        }

        if (false == sf_init_fast_sweep (tau,
                                         n2, n1,
                                         o2, o1,
                                         d2, d1,
                                         sy, sz,fac))
        mexErrMsgTxt("Incorrect shot location");
 
        sf_run_fast_sweep (tau, T0, py0, pz0, 
                           vz, vx, theta,
                           rhs, niter,
                           n2, n1,
                           o2, o1,
                           d2, d1,
                           sy, sz, fac);
 
        for(i=0;i<n12;i++){
 
            if(fac==0) {t[i] = T0[i]+tau[i]; tn1[i] = tn[i]; tn[i] = t[i];}
            else if(fac==1){ t[i] = T0[i]*tau[i]; tn1[i] = tn[i]; tn[i] = t[i];}
            else mexErrMsgTxt("Choose fac=0(Additive) or fac=1(Multiplicative) factorization");
             
        }
 
        sf_derivative_2D(tau,pz,py,d2,d1,n2,n1,2);
 
        sum = 0.;
 
        /* dT/dx and dT/dz computation */
        for(i=0;i<n12;i++){
            if(fac==0){
                py[i] = py[i]+py0[i]; pz[i] = pz[i]+pz0[i];
            }
            else if(fac==1){
                py[i] = T0[i]*py[i] + tau[i]*py0[i]; pz[i] = T0[i]*pz[i] + tau[i]*pz0[i];
            }
            sum = sum + fabs(tn[i]-tn1[i]);
        }
 
       if(info!= 0) mexPrintf("================================== \n");
       if(info!= 0) mexPrintf("L1 norm of update is %f \n",sum/n12);
       if(info!= 0) mexPrintf("================================== \n");
 
        for(i=0;i<n1;i++){
            for(j=0;j<n2;j++){
  
                pydash = cos(theta[j*n1+i])*py[j*n1+i] + sin(theta[j*n1+i])*pz[j*n1+i];
                pzdash = cos(theta[j*n1+i])*pz[j*n1+i] - sin(theta[j*n1+i])*py[j*n1+i];
                rhs[j*n1+i] = 1 + 2*((eps[j*n1+i]-del[j*n1+i])/(1+2*eps[j*n1+i]))
                       *vx[j*n1+i]*vx[j*n1+i]*vz[j*n1+i]*vz[j*n1+i]*pydash*pydash*pzdash*pzdash;
 
            }
        }
    }
        
    //Fill output array
    for (i=0;i<G.nx;i++){
        for (j=0;j<G.nz;j++){
            out[i+j*G.nx] = (double)t[j+i*G.nz];
        }
    }
 
    delete[] t; delete[] tau; delete[] vz; delete[] vx; delete[] eps;
    delete[] del; delete[] theta; delete[] T0; delete[] py0;
    delete[] pz0; delete[] py; delete[] pz; delete[] rhs;
    delete[] tn; delete[] tn1;

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