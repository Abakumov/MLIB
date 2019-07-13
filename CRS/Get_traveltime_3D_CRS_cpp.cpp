// 3D i-CRS stacking operator
// Abakumov Ivan
// University of Hamburg
// e-mail: abakumov_ivan@mail.ru
// 3rd of September 2016
// Hamburg
//////////////////////////////////////////////////////
//c++  libraries
#include "mex.h"
#include <iostream>
#include <cstdio>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <cstdlib>
//using
using namespace std;
using std::cerr;
using std::endl;
using std::ofstream;
//constants
#define info 1

/*
 attr[0]=2/v0*cos(az)*sin(dip)
 attr[1]=2/v0*sin(az)*sin(dip)
 
 attr[2]=nip00                      M(1,1)
 attr[3]=nip10                    2*M(1,2)
 attr[4]=nip11                      M(2,2)
 
 attr[5]=n00                        N(1,1)
 attr[6]=n10                      2*N(1,2)
 attr[7]=n11                        N(2,2)
*/
// 
void CRSstack(int nt, double *attr, double *mx, double *my, double *hx, double *hy, double *mxmx, double *mxmy, double *mymy, double *hxhx, double *hxhy, double *hyhy, double *tti)
{
    int tnum; 
    double t0; 
    
    t0 = attr[8]; 
     
    for (tnum=0;tnum<nt;tnum++)
    {
    /////////////////////////////////////////////////////    
        double t2=t0+attr[0]*mx[tnum]+attr[1]*my[tnum];
        t2*=t2;
        t2+=2*t0*(attr[2]*hxhx[tnum]+2*attr[3]*hxhy[tnum]+attr[4]*hyhy[tnum]);
        t2+=2*t0*(attr[5]*mxmx[tnum]+2*attr[6]*mxmy[tnum]+attr[7]*mymy[tnum]);
        //t2>0?sqrt(t2):-1.0;
        tti[tnum] = sqrt(t2);
    //////////////////////////////////////////////////////       
    }
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*---- Gateway routine-----*/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
//declare variables
const mwSize *dims;
mwSize ndim;
int nt;
double *tti;
double *attr, *mx, *my, *hx, *hy;
double *mxmx, *mxmy, *mymy, *hxhx, *hxhy, *hyhy;

//check input
if (nrhs!=11) mexErrMsgTxt("Three input arguments are required.");
if (nlhs!=1) mexErrMsgTxt("One output argument is required.");
//associate input pointers
attr = (double *)mxGetPr(prhs[0]);
mx   = (double *)mxGetPr(prhs[1]);
my   = (double *)mxGetPr(prhs[2]);
hx   = (double *)mxGetPr(prhs[3]);
hy   = (double *)mxGetPr(prhs[4]);
mxmx = (double *)mxGetPr(prhs[5]);
mxmy = (double *)mxGetPr(prhs[6]);
mymy = (double *)mxGetPr(prhs[7]);
hxhx = (double *)mxGetPr(prhs[8]);
hxhy = (double *)mxGetPr(prhs[9]);
hyhy = (double *)mxGetPr(prhs[10]);
//check if dimensions of velmod are correct
ndim = mxGetNumberOfDimensions(prhs[1]);
dims = mxGetDimensions(prhs[1]);
if (ndim!=2) mexErrMsgTxt("Number of dimensions of velmod array is not correct.");
if (info!= 0) mexPrintf("dims: dims[0]=%d, dims[1]=%d\n",dims[0],dims[1]);
nt=dims[1];
//create array for output
plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
//associate output pointer
tti = (double *)mxGetPr(plhs[0]);
//FSM(size, G, shot, velmod, tti);
CRSstack(nt,attr,mx,my,hx,hy,mxmx,mxmy,mymy,hxhx,hxhy,hyhy,tti);
return;
}
