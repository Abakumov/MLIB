/*
--- Elastic 2D forward modelling code based on pseudo-Fourier method ---
Parameters:
- For standard modelling with standard source, 
    Source type:
    stype  0: pressure source,  1: horizontal force,  2: vertical force,  3: shear
    Wavelet type: 
    wtype  1: Gaussian 2: 1st derivative if Gaussian ,3: Ricker(2nd derivivative of Gaussian), 4: Seismic Ricker
- For modelling with boundary condition as source, 
    Boundary condition (bc) type:
    stype  0: pressure bc, 1: Ux bc, 2: Uz bc, 3: Ux & Uz bc;
------------------------------------------------------------------------
Authors:
Denis Anikiev (denis.anikiev@zmaw.de), Ekkehart Tessmer
Saint Petersburg State University / University of Hamburg (2013-2014)
------------------------------------------------------------------------
Version history:
 v1.0 (05.05.2013): 
 *basic version
 v1.1 (20.05.2013): 
 *damping weights are provided
 v1.2 (22.05.2013): 
 *checks stability criteria
 v2.0 (11.07.2013): 
 *global restructurization
 *parallel version for many shots / boundary conditions
 *moment tensor source type is omitted
 *static variables eliminated (except pfafft.cpp)
 v2.1 (11.07.2013): 
 *several functions are put inside el2d.cpp
 v2.2 (12.07.2013):
 *output of displacement derivatives
 *renamed functions
 *difx, difz, rc are inside
 *datatypes.h unnecessary
 v2.3 (15.07.2013):
 *improved input
 *receivers could be in arbitrary grid positions
 *snaps are saved for local target zone
 *dtype renamed to float, itype substituted by int
 *3 output modes
 v2.4 (16.12.2013)
 *special input case - float ux snapshot
 v2.5 (06.03.2014)
 *removed bug with wavelet scaling
 *downgrade: simplified output - no separation to sections and snapshots, no output modes
 *corrected formulation of pressure
 v2.6 (12.03.2014)
 *seismogram section output revised - ilx array in a form [irx,irz] and ilz=[] will result in section output
 *various number of input arrays is possible (memory is not allocated for unnecessary fields)
 *removed bug in automatic grid filling
 v2.7 (04.04.2014) 
 *improvement of the speed due to precomputation of coefficients
 *optimization of difx and difz routines
 v3.0 (05.04.2014)
 *implementation of FFTW library instead of CWP pfafft
 *FFTW is faster
 *non-ANSI headers and CWP headers are unnecessary
 v3.1 (07.04.2014)
 *improvement of absorbing boundary
 *verification of data class of input arguments: CheckClassName function
 v4.0 (14.04.2014)
 *correction of bug that was related to wrong sourcve and wrong order of field output
 *rearrangement of functions
 *improved absorbing boundary 
 *changed order of a source type
 *changed output order
*/
#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <omp.h>
#include <fftw3.h>

/*---- Definitions ----*/
using std::vector;
using std::max;
using std::min;
using std::strcmp;
typedef vector<int> ivector;
typedef vector<float> fvector;
typedef vector<fvector> fpair;

/*---- Global constants ------*/
#ifndef info//information
#define info 1 
#endif

#ifndef fto//fourth order scheme
#define fto 1
#endif

#ifndef BIG//big value to track stability
#define BIG 1.e10
#endif

/*---- Structures ----*/
struct pmtstruct {//general parameters
  int nx,nz,nt,stype,atype,wtype,nxl,nzl,nsht,nb,ngob,bc,grid,nout;
  float dx,dz,dt,fpeak,alpha,vmax,vmin;
};
struct dftstruct {//parameters for FFTW
  fftwf_plan m_dx, p_dx, m_dz, p_dz;
  fftwf_complex *in_dx;
  fftwf_complex *in_dz;
};

/*---- Auxiliary routines ----*/
/* Check data class name */
void CheckClassName(const mxArray* prhs[], int n, const char *class_name)
{
    const char *arg_class_name = mxGetClassName(prhs[n]);
    if (strcmp(arg_class_name,class_name)) {
        mexPrintf("Wrong class name for input argument %d: %s\n", n+1, arg_class_name);
        mexPrintf("Change class to: %s\n", class_name);
        mexErrMsgTxt(" ");
    }
}
/* Rounding */
int iround(float num) 
{
    return (num > 0.0) ? (int)floor(num + 0.5) : (int)ceil(num - 0.5);
}
float dround(float num) 
{
    return (num > 0.0) ? (float)floor(num + 0.5) : (float)ceil(num - 0.5);
}

/* Damping coefficients for absorbing boundaries */
void gb(fvector &gob, int ngob, float alpha)
{
    float pi = 4.*atan(1.);
    for (int i=0; i<ngob; i++) {        
        // cosine:
        gob[i] = 1.0 - alpha*(1.0 + cos((float)i/(float(ngob))*pi))*0.5;
    }
}

/* Wave numbers */
void rk(fvector &ak, int n, float d, int ind)
{
    int i,isign,n2;
    float pi, dn, c;
    pi = 4.*atan(1.);
    isign = (int)(-pow((float)(-1.),(float)ind));
    n2 = n/2;
    dn = 2.*pi/((float)n*d);
    for (i=0; i<n; i++) {
        c=dn*i;
        if (i > n2) c = -dn*(n-i);
        ak[i] = pow(c,(float)ind)*isign/n;
    }
}

/* Spatial derivative in x */
void difx(const fpair &vin, fpair &vout, dftstruct dft, const fvector &rkx, float dx, int iadd)
{
    int nx   = vin[0].size();
    int nz   = vin.size();
    int ieo = nz % 2;    
    int num = ((nz + 1) / 2) * 2;
    float *a = (float *)dft.in_dx;
    
    /* Load data from vin */
    if (ieo) {
        for (int k=0; k<num-2; k+=2) {
            for (int i=0; i<nx; i++) {
                a[k*nx+2*i]   = vin[k][i];
                a[k*nx+2*i+1] = vin[k+1][i];
            }
        }
        for (int i=0; i<nx; i++) {
            a[nx*(nz-1)+2*i]   = vin[nz-1][i];
            a[nx*(nz-1)+2*i+1] = 0.;
        }
    } else {
        for (int k=0; k<num; k+=2) {
            for (int i=0; i<nx; i++) {
                a[k*nx+2*i]   = vin[k][i];
                a[k*nx+2*i+1] = vin[k+1][i];
            }
        }
    }
    
    /* Perform derivative */
    fftwf_execute(dft.m_dx);
    
    for (int k=0; k<num; k+=2) {
        for (int i=0; i<nx; i++) {
            float tmp    =  a[k*nx+2*i]   * rkx[i];
            a[k*nx+2*i]   = -a[k*nx+2*i+1] * rkx[i];
            a[k*nx+2*i+1] =  tmp;
        }
    }
    
    fftwf_execute(dft.p_dx);
    
    if (iadd == 0) {        
        /* Set data into vout */
        if (ieo) {
            for (int k=0; k<num-2; k+=2) {
                for (int i=0; i<nx; i++) {
                    vout[k][i]   = a[k*nx+2*i];
                    vout[k+1][i] = a[k*nx+2*i+1];
                }
            }
            for (int i=0; i<nx; i++) vout[nz-1][i] = a[nx*(nz-1)+2*i];
        } else {
            for (int k=0; k<num; k+=2) {
                for (int i=0; i<nx; i++) {
                    vout[k][i]   = a[k*nx+2*i];
                    vout[k+1][i] = a[k*nx+2*i+1];
                }
            }
        }
    } else if (iadd == 1){        
        /* Add data into vout */
        if (ieo) {
            for (int k=0; k<num-2; k+=2) {
                for (int i=0; i<nx; i++) {
                    vout[k][i]   += a[k*nx+2*i];
                    vout[k+1][i] += a[k*nx+2*i+1];
                }
            }
            for (int i=0; i<nx; i++) vout[nz-1][i] += a[nx*(nz-1)+2*i];
        } else {
            for (int k=0; k<num; k+=2) {
                for (int i=0; i<nx; i++) {
                    vout[k][i]   += a[k*nx+2*i];
                    vout[k+1][i] += a[k*nx+2*i+1];
                }
            }
        }    
    }
    return;
}

/* Spatial derivative in z */
void difz(const fpair &vin, fpair &vout, dftstruct dft, const fvector &rkz, float dz, int iadd)
{
    int nx   = vin[0].size();
    int nz   = vin.size();
    int ieo = nx % 2;    
    int num = ((nx + 1) / 2) * 2;
    float *a = (float *)dft.in_dz;
    
    /* Load data from vin */
    for (int k=0; k<nz; k++) {
        for (int i=0; i<nx; i++) a[k*num+i] = vin[k][i];
    }
    if (ieo) {
        for (int k=0; k<nz; k++) a[k*num+nx] = 0.;
    }
    
    /* Perform derivative */
    fftwf_execute(dft.m_dz);
    
    for (int k=0; k<nz; k++) {
        for (int i=0; i<num; i+=2) {
            float tmp   =  a[k*num+i]   * rkz[k];
            a[k*num+i]   = -a[k*num+i+1] * rkz[k];
            a[k*num+i+1] =  tmp;
        }
    }
    
    fftwf_execute(dft.p_dz);
    
    if (iadd == 0) {        
        /* Set data into vout */
        for (int k=0; k<nz; k++) {
            for (int i=0; i<nx; i++) vout[k][i] = a[k*num+i];
        }
    } else if (iadd == 1) {        
        /* Add data into vout */
        for (int k=0; k<nz; k++) {
            for (int i=0; i<nx; i++) vout[k][i] += a[k*num+i];
        }
    }
    return;
}

/* Apply 2d absorbing boundary conditions */
void ftbc2d(fpair &a, const fvector &gob)
{
    /*
     * special treatment of corners for avoiding double reduction
     * absorbing zones in x- and z-directions have same size
    */
    int nx   = a[0].size();
    int nz   = a.size();    
    int ngob = gob.size();   
    
    /* treatment of sides */
    #if 1
    for (int k=0; k<nz; k++) {
        int ng = ngob-1;
        if (k < ngob) ng = k;
        if (k > nz-ngob) ng=nz-k-1;
        for (int ig=0; ig<ng+1; ig++) {            
            float g = gob[ig];            
            a[k][ig]      *= g;
            a[k][nx-ig-1] *= g;
        }
    }
    #endif
    /* treatment of top and bottom */
    #if 1
    for (int k=0; k<ngob; k++) {
        int ig = k;
        int ia = ig+1;
        int iex = nx - ig -1;        
        for (int ii=ia; ii<iex; ii++) {            
            float g = gob[ig];
            a[k][ii]      *= g;
            a[nz-k-1][ii] *= g;
        }
    }
    #endif
}

/* Construct elastic parameters 1/rho, 2*mu and lambda from rho,vp,vs */
void parameters(int nx, int nz, float rho[], float vp[], float vs[], fpair &lambda, fpair &tmu, fpair &rhoinv)
{    
    for (int k=0; k<nz; k++) {
        for (int i=0; i<nx; i++) {
            int index = i + k*nx;
            rhoinv[k][i] = 1./rho[index];                               // inverse density
            tmu[k][i] = 2.*vs[index]*vs[index]*rho[index];              // double mu
            lambda[k][i] = vp[index]*vp[index]*rho[index] - tmu[k][i];  //lambda
        }
    }
}

/* Return max amplitude */
float maxamp(int nx, int nz, const fpair &ux, const fpair &uz) 
{
    float amx = 0.;
    float amz = 0.;
    for (int iz=0; iz<nz; iz++) {
        for (int ix=0; ix<nx; ix++) {
            amx = (fabs(ux[iz][ix]) > amx ? fabs(ux[iz][ix]) : amx);
            amz = (fabs(uz[iz][ix]) > amz ? fabs(uz[iz][ix]) : amz);
        }
    }
    return max(amx,amz);
}

/* Wavelet function */
float fwave(int wtype, float t, float fpeak)
{
    float pi,pi2,agauss,tcut,agauss2;
    float s,tmp,s2;   
    
    pi      = 4.*atan(1.);
    pi2     = pi*pi;
    agauss  = fpeak;
    agauss2 = agauss*agauss;
    tcut    = 1.5/agauss;
    s   = (t-tcut)*agauss;
    s2  = s*s;
    tmp = 0;
    
    if (fabs(s) < 4.) {
        switch (wtype) {
            case 1: //gaussian
                tmp = exp(-pi2*s2);
                break;
            case 2: //1st gaussian derivative
                tmp = -2*pi2*agauss*s*exp(-pi2*s2);
                break;
            case 3: //2nd gaussian derivative (ricker)
                tmp = 2*pi2*agauss2*(2*pi2*s2-1)*exp(-pi2*s2);
                break;
            case 4: //seismic ricker wavelet
                tmp = (1-2*pi2*s2)*exp(-pi2*s2);
                break;
            default:
                tmp = 0;
                break;
        }
    }
    return tmp;
}

/* Accellerations */
void acceleration(float dx, float dz, dftstruct dft, const fvector &rkx, const fvector &rkz, fpair &a1, fpair &a2, const fpair &a3)
{
    int nx = a1[0].size();
    int nz = a1.size();
    difx(a1,a1,dft,rkx,dx,0); //d(sxx)/dx
    difz(a3,a1,dft,rkz,dz,1); //d(sxx)/dx + d(sxx)/dz
    difz(a2,a2,dft,rkz,dz,0); //d(szz)/dz
    difx(a3,a2,dft,rkx,dx,1); //d(szz)/dz + d(szz)/dx
    return;
}
/*---- Master routine  -------*/
void el2d(pmtstruct pmt, float rho[], float vp[], float vs[], int sxi[], int szi[], int lxi[], int lzi[], float wavelet[], float damp[], float pb[], float uxb[], float uzb[], float *psn, float *uxsn, float *uzsn, float *uxxsn, float *uzzsn, float *uxzsn, float *uzxsn)
{
    /* --- Parallel loop variable --- */
    int isht;
    
    /* --- Parameters --- */
    int nout     = pmt.nout;
    int grid     = pmt.grid;
    int bc       = pmt.bc;
    int stype    = pmt.stype;
    int atype    = pmt.atype;
    int wtype    = pmt.wtype;
    int nx       = pmt.nx;
    int nz       = pmt.nz;
    int nt       = pmt.nt;
    int nsht     = pmt.nsht;
    int nb       = pmt.nb;
    int nxl      = pmt.nxl;
    int nzl      = pmt.nzl;
    int ngob     = pmt.ngob;
    float dx    = pmt.dx;
    float dz    = pmt.dz;
    float dt    = pmt.dt;
    float fpeak = pmt.fpeak;
    float alpha = pmt.alpha;
    float vmax  = pmt.vmax;
    float dt2   = dt*dt;
    float dt4t  = dt2*dt2/(float)12;
    
    /* --- Arrays --- */
    fvector gob(ngob,0.);
    fpair lambda(nz,fvector(nx,0.));
    fpair tmu(nz,fvector(nx,0.));
    fpair rhoinv(nz,fvector(nx,0.));
    
    /* --- Elastic parameters --- */
    parameters(nx,nz,rho,vp,vs,lambda,tmu,rhoinv);
    
    /* --- Get damping coefficients for absorbing boundaries --- */
    switch (atype) {
        case 1:
            gb(gob,ngob,alpha);
            break;
        case 2:
            for (int k=0; k<ngob; k++) gob[k] = damp[k];            
            break;
    }
    
    /* --- Precompute coefficients --- */
    fvector rkx(nx,0), rkz(nz,0);
    rk(rkx,nx,dx,1);
    rk(rkz,nz,dz,1);
    
    /* --- Parallel section --- */
    #pragma omp parallel for private(isht)
    for (isht=0; isht<nsht; isht++) {
        
        int isx,isz;
        float fw = 0;
        fpair uxp(nz,fvector(nx,0.));
        fpair uzp(nz,fvector(nx,0.));
        fpair ux (nz,fvector(nx,0.));
        fpair uz (nz,fvector(nx,0.));
        fpair p  (nz,fvector(nx,0.));
        fpair uxx(nz,fvector(nx,0.));
        fpair uzz(nz,fvector(nx,0.));
        fpair uxz(nz,fvector(nx,0.));
        fpair uzx(nz,fvector(nx,0.));
        fpair a1 (nz,fvector(nx,0.));
        fpair a2 (nz,fvector(nx,0.));
        fpair a3 (nz,fvector(nx,0.));
        
        /* --- DFT plan variables --- */    
        dftstruct dft;
        int n[1], inembed[1], onembed[1];
        int num, idist, istride, odist, ostride;
        
        /* --- Precompute DFT plans for difx --- */
        num = ((nz + 1) / 2) * 2;
        n[0] = inembed[0] = onembed[0] = nx;
        istride = ostride = 1;
        idist = odist = nx;
        #pragma omp critical
        {
            dft.in_dx = (fftwf_complex *) fftwf_malloc(nx*num/2*sizeof(fftwf_complex));
            dft.m_dx = fftwf_plan_many_dft(1, n, num/2,
                    dft.in_dx, inembed, istride, idist,
                    dft.in_dx, onembed, ostride, odist,
                    FFTW_FORWARD, FFTW_MEASURE);
            dft.p_dx = fftwf_plan_many_dft(1, n, num/2,
                    dft.in_dx, inembed, istride, idist,
                    dft.in_dx, onembed, ostride, odist,
                    FFTW_BACKWARD, FFTW_MEASURE);
        }
        /* --- Precompute DFT plans for difz --- */
        num = ((nx + 1) / 2) * 2;
        n[0] = inembed[0] = onembed[0] = nz;
        istride = ostride = num/2;
        idist = odist = 1;
        #pragma omp critical
        {
            dft.in_dz = (fftwf_complex *) fftwf_malloc(nz*num/2*sizeof(fftwf_complex));
            dft.m_dz = fftwf_plan_many_dft(1, n, num/2,
                    dft.in_dz, inembed, istride, idist,
                    dft.in_dz, onembed, ostride, odist,
                    FFTW_FORWARD, FFTW_MEASURE);
            dft.p_dz = fftwf_plan_many_dft(1, n, num/2,
                    dft.in_dz, inembed, istride, idist,
                    dft.in_dz, onembed, ostride, odist,
                    FFTW_BACKWARD, FFTW_MEASURE);
        }
                        
        /* --- Loop over time steps --- */
        for (int it=0; it<nt; it++) {
            
            /* --- Strains --- */
            difx(ux,uxx,dft,rkx,dx,0); // d(ux)/dx = exx
            difz(uz,uzz,dft,rkz,dz,0); // d(uz)/dz = ezz
            difz(ux,uxz,dft,rkz,dz,0); // d(ux)/dz
            difx(uz,uzx,dft,rkx,dx,0); // d(uz)/dx
            
            /* --- Stress and pressure --- */
            for (int iz=0; iz<nz; iz++) {
                for (int ix=0; ix<nx; ix++) {
                    float u1 = uxx[iz][ix] + uzz[iz][ix];                     // exx + ezz
                    float u2 = uxz[iz][ix] + uzx[iz][ix];                     // d(uz)/dx + d(ux)/dz = 2*exz = 2*ezx
                    a1[iz][ix] = lambda[iz][ix]*u1 + tmu[iz][ix]*uxx[iz][ix]; // sxx = lambda*(exx+ezz) + 2*mu*exx
                    a2[iz][ix] = lambda[iz][ix]*u1 + tmu[iz][ix]*uzz[iz][ix]; // szz = lambda*(exx+ezz) + 2*mu*ezz
                    a3[iz][ix] = tmu[iz][ix]*u2*0.5;                          // 2*mu*exz
                    p[iz][ix]  = -(lambda[iz][ix] + tmu[iz][ix])*u1;          // pressure p = -(lambda+2*mu)*(exx+ezz)
                }
            }

            /* --- Output wavefields --- */            
            #pragma omp critical
            {
                if (grid) {
                    for (int izl=0; izl<nzl; izl++) {
                        int iz = lzi[izl];
                        for (int ixl=0; ixl<nxl; ixl++) {
                            int ix = lxi[ixl];
                            int index  = ixl + izl*nxl + it*nxl*nzl + isht*nt*nxl*nzl;
                            if (nout==7) {
                                psn  [index] = p  [iz][ix];
                                uxsn [index] = ux [iz][ix];
                                uzsn [index] = uz [iz][ix];
                                uxxsn[index] = uxx[iz][ix];
                                uzzsn[index] = uzz[iz][ix];
                                uxzsn[index] = uxz[iz][ix];
                                uzxsn[index] = uzx[iz][ix];
                            } else {
                                if (nout>0) psn  [index] = p  [iz][ix];
                                if (nout>1) uxsn [index] = ux [iz][ix];
                                if (nout>2) uzsn [index] = uz [iz][ix];
                                if (nout>3) uxxsn[index] = uxx[iz][ix];
                                if (nout>4) uzzsn[index] = uzz[iz][ix];
                                if (nout>5) uxzsn[index] = uxz[iz][ix];
                                if (nout>6) uzxsn[index] = uzx[iz][ix];
                            }
                        }
                    }
                } else {
                    for (int irec=0; irec<nxl; irec++) {
                        int ix = lxi[irec];
                        int iz = lzi[irec];
                        int index  = irec + it*nxl + isht*nt*nxl;
                        if (nout==7) {
                            psn  [index] = p  [iz][ix];
                            uxsn [index] = ux [iz][ix];
                            uzsn [index] = uz [iz][ix];
                            uxxsn[index] = uxx[iz][ix];
                            uzzsn[index] = uzz[iz][ix];
                            uxzsn[index] = uxz[iz][ix];
                            uzxsn[index] = uzx[iz][ix];
                        } else {
                            if (nout>0) psn  [index] = p  [iz][ix];
                            if (nout>1) uxsn [index] = ux [iz][ix];
                            if (nout>2) uzsn [index] = uz [iz][ix];
                            if (nout>3) uxxsn[index] = uxx[iz][ix];
                            if (nout>4) uzzsn[index] = uzz[iz][ix];
                            if (nout>5) uxzsn[index] = uxz[iz][ix];
                            if (nout>6) uzxsn[index] = uzx[iz][ix];
                        }
                    }
                }
            }
            
            /* Add source term */
            switch (bc) {
                case 0: // point source
                    isx = sxi[isht];
                    isz = szi[isht];
                    switch (wtype) {
                        case 0://user specified wavelet
                            fw = wavelet[it]/(dx*dz);    
                            break;
                        case 1://construct wavelet
                            float time = (float)it*dt;
                            fw = fwave(wtype,time,fpeak)/(dx*dz);
                            break;
                    }
                    /* --- Apply source term and compute accelerations --- */
                    switch (stype) {
                        case 0: // pressure source
                            a1[isz][isx] -= fw;
                            a2[isz][isx] -= fw;   
                            acceleration(dx,dz,dft,rkx,rkz,a1,a2,a3);
                            break;
                        case 1: // horizontal force source
                            acceleration(dx,dz,dft,rkx,rkz,a1,a2,a3);
                            a1[isz][isx] += fw;
                            break;
                        case 2: // vertical force source
                            acceleration(dx,dz,dft,rkx,rkz,a1,a2,a3);
                            a2[isz][isx] += fw;
                            break;
                        case 3: // shear source
                            a3[isz][isx] += fw;
                            acceleration(dx,dz,dft,rkx,rkz,a1,a2,a3);
                            break;
                    }
                    break;
                case 1: // boundary conditions
                    /* --- Apply boundary condition --- */
                    switch (stype) {
                        case 0: // pressure as boundary condition
                            for (int ib=0; ib<nb; ib++) {
                                int ibx = sxi[ib];
                                int ibz = szi[ib];
                                int index = ib + it*nb + isht*nt*nb;
                                a1[ibz][ibx] -= pb[index];
                                a2[ibz][ibx] -= pb[index];
                            }
                            acceleration(dx,dz,dft,rkx,rkz,a1,a2,a3);
                            break;
                        case 1: // horizontal displacement as boundary condition
                            acceleration(dx,dz,dft,rkx,rkz,a1,a2,a3);
                            for (int ib=0; ib<nb; ib++) {
                                int ibx = sxi[ib];
                                int ibz = szi[ib];
                                int index = ib + it*nb + isht*nt*nb;                                
                                ux[ibz][ibx] += uxb[index];
                            }
                            break;
                        case 2: // vertical displacement as boundary condition
                            acceleration(dx,dz,dft,rkx,rkz,a1,a2,a3);
                            for (int ib=0; ib<nb; ib++) {
                                int ibx = sxi[ib];
                                int ibz = szi[ib];
                                int index = ib + it*nb + isht*nt*nb;
                                uz[ibz][ibx] += uzb[index];
                            }
                            break;
                        case 3: // both horizontal and vertical displacement as boundary condition
                            acceleration(dx,dz,dft,rkx,rkz,a1,a2,a3);
                            for (int ib=0; ib<nb; ib++) {
                                int ibx = sxi[ib];
                                int ibz = szi[ib];
                                int index = ib + it*nb + isht*nt*nb;
                                ux[ibz][ibx] += uxb[index];
                                uz[ibz][ibx] += uzb[index];
                            }
                            break;
                    }
                    break;
            }

            /* --- Multiply by the inverse of density --- */
            for (int iz=0; iz<nz; iz++) {
                for (int ix=0; ix<nx; ix++) {
                    a1[iz][ix] *= rhoinv[iz][ix];
                    a2[iz][ix] *= rhoinv[iz][ix];
                }
            }

            /* --- Apply absorbing boundaries --- */
            if (atype > 0) {
                ftbc2d(a1,gob);
                ftbc2d(a2,gob);
            }

            /* --- Time integration --- */
            for (int iz=0; iz<nz; iz++) {
                for (int ix=0; ix<nx; ix++) {
                    uxp[iz][ix] = 2.*ux[iz][ix] - uxp[iz][ix] + dt2*a1[iz][ix];
                    uzp[iz][ix] = 2.*uz[iz][ix] - uzp[iz][ix] + dt2*a2[iz][ix];
                }
            }
            
            /* --- Swap arrays --- */
            for (int iz=0; iz<nz; iz++) {
                for (int ix=0; ix<nx; ix++) {
                    float u1   = uxp[iz][ix];
                    float u2   = uzp[iz][ix];
                    uxp[iz][ix] = ux[iz][ix];
                    uzp[iz][ix] = uz[iz][ix];
                    ux[iz][ix]  = u1;
                    uz[iz][ix]  = u2;
                }
            }

            /* --- Fourth order extension ----*/
            if (fto) {
                difz(a1,a3,dft,rkz,dz,0);
                difx(a2,a3,dft,rkx,dx,1);
                difx(a1,a1,dft,rkx,dx,0);
                difz(a2,a2,dft,rkz,dz,0);
                
                for (int k=0; k<nz; k++) {
                    for (int i=0; i<nx; i++) {
                        float theta = a1[k][i] + a2[k][i];
                        a1[k][i] = lambda[k][i]*theta + tmu[k][i]*a1[k][i];
                        a2[k][i] = lambda[k][i]*theta + tmu[k][i]*a2[k][i];
                        a3[k][i] = tmu[k][i]*a3[k][i]*0.5;
                    }
                }
                
                difx(a1,a1,dft,rkx,dx,0);
                difz(a3,a1,dft,rkz,dz,1);                
                difz(a2,a2,dft,rkz,dz,0);
                difx(a3,a2,dft,rkx,dx,1);
                
                if (atype > 0) {
                    ftbc2d(a1,gob);
                    ftbc2d(a2,gob);                    
                }
                
                for (int k=0; k<nz; k++) {
                    for (int i=0; i<nx; i++) {
                        ux[k][i] += a1[k][i] * rhoinv[k][i] * dt4t;
                        uz[k][i] += a2[k][i] * rhoinv[k][i] * dt4t;
                    }
                }
            }
            
            /* --- apply absorbing boundaries --- */
            if (atype > 0) {
                ftbc2d(ux ,gob);
                ftbc2d(uz ,gob);
                ftbc2d(uxp,gob);
                ftbc2d(uzp,gob);
            }
            
            /* --- check stability --- */
            if (it%50 == 0) {
                float amax = maxamp(nx,nz,ux,uz);
                if (amax > BIG) {
                    #pragma omp critical
                    {                        
                        printf("Unstable run: shot %4d; step %4d; amax=%e;\n",isht,it,amax);
                    }
                } else if (info) {
                    #pragma omp critical
                    {
                        printf("shot %4d; step %4d; amax=%e;\n",isht,it,amax);
                    }
                }
            }
        } // end of loop over time steps
        
        #pragma omp critical
        { 
            /* --- Clean up and free memory --- */
            fftwf_destroy_plan(dft.m_dx);
            fftwf_destroy_plan(dft.p_dx);
            fftwf_destroy_plan(dft.m_dz);
            fftwf_destroy_plan(dft.p_dz);
            fftwf_free(dft.in_dx);
            fftwf_free(dft.in_dz);
        }        
    }
    /* End of parallel for loop */
    
    fftwf_cleanup();    
    return;
}

/*----------------------------------------------------------------------*/
/*---- Gateway routine -----*/
/*----------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    pmtstruct pmt;
    mxClassID category;
    int nelem,nsx,nsz,nwav;
    float pi = 4.*atan(1.);
    float *pb = NULL, *uxb = NULL, *uzb = NULL;    //input boundary conditions    
    float *psn = NULL, *uxsn = NULL, *uzsn = NULL; //output snapshots    
    float *uxxsn, *uzzsn, *uxzsn, *uzxsn;          //
    float *rho = NULL, *vp = NULL, *vs = NULL;     //elastic parameters
    float *wavelet = NULL, *damp = NULL;           //wavelet and damping coefficients arrays
    int *sxi = NULL, *szi = NULL;                   //shot index arrays
    int *lxi = NULL, *lzi = NULL;                   //target local zone index arrays
    int *ri = NULL;                                 //receiver index array
    ivector xi, zi;                                 //extra arrays
    
    /* Check for proper number of input and output arguments */
    if (nrhs>17 || nrhs<15) mexErrMsgTxt("Incorrect number of input arguments!");
    
    if (nlhs<1 || nlhs>7) { 
        mexErrMsgTxt("Incorrect number of output arguments!");
    } else {
        pmt.nout = nlhs;
    }
    
    /* Boundary condition flag (controls if source is of a standard type or boundary condition) */
    pmt.bc = iround(mxGetScalar(prhs[0]));
    
    /* Model */
    for (int i=1;i<6;i++) CheckClassName(prhs,i,"single");
    rho = (float*)mxGetData(prhs[1]);
    vp  = (float*)mxGetData(prhs[2]);
    vs  = (float*)mxGetData(prhs[3]);
    pmt.nx = (int)mxGetM(prhs[1]);
    pmt.nz = (int)mxGetN(prhs[1]);
    pmt.dx = (float)mxGetScalar(prhs[4]);
    pmt.dz = (float)mxGetScalar(prhs[5]);
    
    /* Absorbing boundary parameters */
    if (mxIsEmpty(prhs[6])) {
        pmt.atype = 0;
        pmt.ngob  = 1;
        pmt.alpha = 0;        
    } else {
        CheckClassName(prhs,6,"single");
        damp = (float*)mxGetData(prhs[6]);
        nelem = (int)mxGetNumberOfElements(prhs[6]);
        if (nelem == 1) {
            pmt.atype = 1;
            pmt.ngob  = iround(damp[0]);
            pmt.alpha = 0.075; //default value
        } else if (nelem == 2) {
            pmt.atype = 1;
            pmt.ngob  = iround(damp[0]);
            pmt.alpha = damp[1];
        } else {
            pmt.atype = 2;
            pmt.ngob  = nelem;
            pmt.alpha = 0;            
        }
    }
    
    /* Read coordinates of local grid for output */
    if (mxIsEmpty(prhs[7])) {
        pmt.grid = 1;
        pmt.nxl = pmt.nx;
        pmt.nzl = pmt.nz;
        xi.assign(pmt.nx,0); 
        zi.assign(pmt.nz,0);
        for (int ix=0; ix<pmt.nx; ix++) xi[ix]=ix;
        for (int iz=0; iz<pmt.nz; iz++) zi[iz]=iz;
        lxi = &xi[0];
        lzi = &zi[0];
    } else if (mxIsEmpty(prhs[8])) {
        if (mxGetM(prhs[7])!=2) mexErrMsgTxt("Receiver coordinates are in wrong format!");
        pmt.grid = 0;
        int nrec = mxGetN(prhs[7]);
        pmt.nxl = nrec;
        pmt.nzl = 1;
        CheckClassName(prhs,7,"int32");
        ri = (int*)mxGetData(prhs[7]);        
        xi.assign(nrec,0);
        zi.assign(nrec,0);
        for (int irec=0; irec<nrec; irec++) {
            xi[irec] = ri[2*irec];
            zi[irec] = ri[2*irec+1];
        }
        lxi = &xi[0];
        lzi = &zi[0];      
    } else {
        pmt.grid = 1;
        CheckClassName(prhs,7,"int32");
        CheckClassName(prhs,8,"int32");        
        pmt.nxl = (int)mxGetNumberOfElements(prhs[7]);
        pmt.nzl = (int)mxGetNumberOfElements(prhs[8]);
        lxi = (int*)mxGetData(prhs[7]);
        lzi = (int*)mxGetData(prhs[8]);
    }
    /* Source type (or boundary condition type) */
    CheckClassName(prhs,9,"int32");
    pmt.stype = (int)mxGetScalar(prhs[9]);
    
    /* Coordinates of shots or nodes where boundary condition is applied */
    CheckClassName(prhs,10,"int32");    
    CheckClassName(prhs,11,"int32");
    nsx = (int)mxGetNumberOfElements(prhs[10]);
    nsz = (int)mxGetNumberOfElements(prhs[11]);
    if (nsx==nsz) {
        sxi = (int*)mxGetData(prhs[10]);
        szi = (int*)mxGetData(prhs[11]);        
    } else {
        mexErrMsgTxt("Numbers of shot x and z coordinates are not consistent!");
    }        
    
    /* Peak frequency of source */
    CheckClassName(prhs,12,"single");    
    pmt.fpeak = (float)mxGetScalar(prhs[12]);
    
    /* Time step */
    CheckClassName(prhs,13,"single");
    pmt.dt = (float)mxGetScalar(prhs[13]);
    
    /* Number of time steps */
    CheckClassName(prhs,14,"int32");    
    pmt.nt = (int)mxGetScalar(prhs[14]);
    
    /* Other parameters depend on source representation type */
    switch (pmt.bc) {
        case 0: /* Standard point source, no boundary condition*/
            pmt.nsht = nsx; //number of shots
            pmt.nb   = 0;            
            nwav = (int)mxGetNumberOfElements(prhs[15]);
            if (nwav==1) {                
                CheckClassName(prhs,15,"int32");
                pmt.wtype = (int)mxGetScalar(prhs[15]);
            } else if (pmt.nt==nwav) {
                pmt.wtype = 0;                
                CheckClassName(prhs,15,"single");
                wavelet   = (float*)mxGetData(prhs[15]);
            } else {
                mexErrMsgTxt("Number of time steps is not consistent with wavelet length!");
            }
            break;
        case 1: /* Boundary condition */            
            pmt.wtype  = 0;            
            pmt.nb     = nsx; //number of nodes for boundary condition
            CheckClassName(prhs,15,"single");
            nelem  = (int)mxGetNumberOfElements(prhs[15]);
            pmt.nsht = nelem/pmt.nb/pmt.nt; //number of recorded boundary conditions
            switch (pmt.stype) {
                case 0: //pressure boundary condition
                    pb  = (float*)mxGetData(prhs[15]);
                    break;
                case 1: //horizontal displacement boundary condition
                    uxb = (float*)mxGetData(prhs[15]);
                    break;
                case 2: //vertical displacement boundary condition
                    uzb = (float*)mxGetData(prhs[15]);
                    break;
                case 3: //horizontal and vertical boundary conditions
                    uxb = (float*)mxGetData(prhs[15]);
                    if (nrhs>15) {                                              // IVAN >18 in initial code
                        if (!mxIsEmpty(prhs[16])) {                             // Do we really need ! here???
                            CheckClassName(prhs,16,"single");
                            uzb = (float*)mxGetData(prhs[16]);
                        } else {
                            mexErrMsgTxt("Missing uz boundary condition!");
                        }
                    } else {
                        mexErrMsgTxt("Missing uz boundary condition!");
                    }
                    break;
                default:
                    mexErrMsgTxt("Unknown boundary condition type!");
                    break;
            }
            break;
        default:
            mexErrMsgTxt("Wrong modelling type!\n");
            break;
    }
    
    /* Set dimensions of output arrays */
    mwSize ndimg = 4, dimg[4]; // grid output
    mwSize ndims = 3, dims[3]; // section output
    dimg[0] = pmt.nxl; dimg[1] = pmt.nzl; dimg[2] = pmt.nt; dimg[3] = pmt.nsht;
    dims[0] = pmt.nxl; dims[1] = pmt.nt;  dims[2] = pmt.nsht;
    /* Create output arrays */
    if (pmt.grid) {
        if (info) mexPrintf("Dimensions for output arrays: %d x %d x %d x %d\n",dimg[0],dimg[1],dimg[2],dimg[3]);
        for (int ilhs=0; ilhs<nlhs; ilhs++) {
            plhs[ilhs]  = mxCreateNumericArray(ndimg, dimg, mxSINGLE_CLASS, mxREAL);
        }    
    } else {
        if (info) mexPrintf("Dimensions for output arrays: %d x %d x %d\n",dims[0],dims[1],dims[2]);
        for (int ilhs=0; ilhs<nlhs; ilhs++) {
            plhs[ilhs]  = mxCreateNumericArray(ndims, dims, mxSINGLE_CLASS, mxREAL);
        }
    }
    if (nlhs>0) psn   = (float*)mxGetData(plhs[0]);
    if (nlhs>1) uxsn  = (float*)mxGetData(plhs[1]);
    if (nlhs>2) uzsn  = (float*)mxGetData(plhs[2]);
    if (nlhs>3) uxxsn = (float*)mxGetData(plhs[3]);
    if (nlhs>4) uzzsn = (float*)mxGetData(plhs[4]);
    if (nlhs>5) uxzsn = (float*)mxGetData(plhs[5]);
    if (nlhs>6) uzxsn = (float*)mxGetData(plhs[6]);    
    
    /* Check stability criteria */
    pmt.vmin = vs[0];
    pmt.vmax = vp[0];
    for (int i=0; i<pmt.nx*pmt.nz; i++) {        
        if (min(vs[i],vp[i]) < pmt.vmin) pmt.vmin = min(vs[i],vp[i]);
        if (max(vs[i],vp[i]) > pmt.vmax) pmt.vmax = max(vs[i],vp[i]);        
    }
    float gamma = (float)2.;
    if (fto) gamma = sqrt((float)12.);
    float crit1 = (float)0.5*(pmt.vmin/pmt.fpeak)/max(pmt.dx,pmt.dz);
    float crit2 = gamma/(pmt.dt*pi*pmt.vmax*sqrt(1./(pmt.dx*pmt.dx)+1./(pmt.dz*pmt.dz)));
    if (crit1 >= 2.) {
        if (info) mexPrintf("\nMinimum number of points per shortest wavelength: %g >= 2.\n",crit1);
    } else {
        mexPrintf("\nWarning: Minimum number of points per shortest wavelength: %g < 2!\n",crit1);
    }
    if (crit2 > 1.) {
        if (info) mexPrintf("\nStability parameter: %g > 1.\n",crit2);
    } else {
        mexPrintf("\nWarning: stability parameter: %g <= 1!\n",crit2);
    }
    
    /* Info about parameters */
    if (info) {
        mexPrintf("Parameters:\n");
        mexPrintf("nout = %d;\n",pmt.nout);
        mexPrintf("nx = %d; nz = %d;\n",pmt.nx,pmt.nz);
        mexPrintf("nxl = %d; nzl = %d;\n",pmt.nxl,pmt.nzl);
        mexPrintf("nsht = %d; nb = %d; ngob=%d;\n", pmt.nsht,pmt.nb,pmt.ngob);
        mexPrintf("nt = %d; dt = %f; fpeak = %f;\n",pmt.nt,pmt.dt,pmt.fpeak);
        mexPrintf("dx = %f; dz = %f;\n",pmt.dx,pmt.dz);
        mexPrintf("bc = %d, stype = %d; wtype = %d; atype=%d;\n",pmt.bc,pmt.stype,pmt.wtype,pmt.atype);
        //if (pmt.wtype==0) mexPrintf("wavelet[0] = %f;\n",wavelet[0]);        
        mexPrintf("sxi[0] = %d; szi[0] = %d;\n",sxi[0],szi[0]);        
        mexPrintf("lxi[0] = %d; lzi[0] = %d;\n",lxi[0],lzi[0]);
    }
    
    /* Call master routine */
    if (info) mexPrintf("Calling master routine from gateway routine...\n");
    el2d(pmt,rho,vp,vs,sxi,szi,lxi,lzi,wavelet,damp,pb,uxb,uzb,psn,uxsn,uzsn,uxxsn,uzzsn,uxzsn,uzxsn);
    if (info) mexPrintf("Master routine terminated successfully.\n");
}