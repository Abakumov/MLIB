// 3D iCRS stacking operator
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
#define info 0

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
void iCRSstack(int nt, double *attr, double *mx, double *my, double *hx, double *hy, double *mxmx, double *mxmy, double *mymy, double *hxhx, double *hxhy, double *hyhy, double *tti)
{
    int tnum; 
    double t0; 
    
    t0 = attr[8]; 
     
    for (tnum=0;tnum<nt;tnum++)
    {
    /////////////////////////////////////////////////////    
        double dxs = mx[tnum]-hx[tnum];
        double dys = my[tnum]-hy[tnum];
        double dxg = mx[tnum]+hx[tnum];
        double dyg = my[tnum]+hy[tnum];
        
        double sNs = attr[5]*dxs*dxs+2*attr[6]*dxs*dys+attr[7]*dys*dys;
        double gNg = attr[5]*dxg*dxg+2*attr[6]*dxg*dyg+attr[7]*dyg*dyg;
        double hMNh = (attr[2]-attr[5])*hxhx[tnum]+2*(attr[3]-attr[6])*hxhy[tnum]+(attr[4]-attr[7])*hyhy[tnum];
        
        double ts0 = t0+attr[0]*dxs+attr[1]*dys;
        double ts1 = ts0*ts0 + 2*t0*sNs;

        double tg0 = t0+attr[0]*dxg+attr[1]*dyg;
        double tg1 = tg0*tg0 + 2*t0*gNg;
        
        double tsg = (sqrt(ts1)+sqrt(tg1))/2;
        
        //  return min(ts1,tg1)>0?0.5*(sqrt(ts1)+sqrt(tg1)):-1.0;
        tti[tnum] = sqrt( tsg*tsg + 2*t0*hMNh);
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
iCRSstack(nt,attr,mx,my,hx,hy,mxmx,mxmy,mymy,hxhx,hxhy,hyhy,tti);
return;
}