// Eikonal Equation solver by Fast Sweeping Method
// 2D version for reflected waves
// Abakumov Ivan
// St. Petersburg University
// e-mail: abakumov_ivan@mail.ru
// 3rd of December 2013
//////////////////////////////////////////////////////
//c++ libraries
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
#define np 3 // Number of points around the shot where vel is "constant". default value 2 or 3

/*---- Auxiliary routines ----*/
double minvalue(double* array, int size)
{
  static double val = array[0];
  for (int i = 1; i < size; i++) {
    val = val <= array[i] ? val : array[i];
  }
  return val;
}

// //Exact traveltime routine for zone around source
double traveltime(double t0, double *x1, double *x2, double v)
{
  double t;
  t = t0 + sqrt(pow(x1[0]-x2[0],2)+pow(x1[1]-x2[1],2))/v;
  return(t);
}

// Routine for evaluation of time using FSM formulae
double xeval(double a, double b, double d)
{
  double x,x1,x2,roots[2];
  roots[0]=min(a,b); roots[1]=max(a,b);
  a = roots[0]; b = roots[1];
  x1 = a + d;
  if (x1 < b) {
    x = x1;
  } else {
    x = (a + b + sqrt(2*pow(d,2) - pow((a-b),2)))/2;
  }
  return(x);
}

// FSM algorithm routine
void fsmalgorithm(int *size, double h, double **vel, double **u)
{
  int i,j,q; // Loop variables
  int xmin,xmax,zmin,zmax; // Index limits
  double a,b,c,d; // Temporary
  xmin=1; xmax=size[0]+1; zmin=1; zmax=size[1]+1;
  for (q=0;q<2;q++) { //2 iterations of FSM
    //1
    for (i=xmin;i<xmax;i++) {
      for (j=zmin;j<zmax;j++) {
        a = min(u[i-1][j], u[i+1][j]);
        b = min(u[i][j-1], u[i][j+1]);
        d = h/vel[i][j];
        u[i][j] = min( u[i][j], xeval(a,b,d) );
      }
    }
    //2
    for (i=xmin;i<xmax;i++) {
      for (j=zmax-1;j>zmin-1;j--) {
        a = min(u[i-1][j], u[i+1][j]);
        b = min(u[i][j-1], u[i][j+1]);
        d = h/vel[i][j];
        u[i][j] = min( u[i][j], xeval(a,b,d) );
      }
    }
    //3
    for (i=xmax-1;i>xmin-1;i--) {
      for (j=zmin;j<zmax;j++) {
        a = min(u[i-1][j], u[i+1][j]);
        b = min(u[i][j-1], u[i][j+1]);
        d = h/vel[i][j];
        u[i][j] = min( u[i][j], xeval(a,b,d) );
      }
    }
    //4
    for (i=xmax-1;i>xmin-1;i--) {
      for (j=zmax-1;j>zmin-1;j--) {
        a = min(u[i-1][j], u[i+1][j]);
        b = min(u[i][j-1], u[i][j+1]);
        d = h/vel[i][j];
        u[i][j] = min( u[i][j], xeval(a,b,d) );
      }
    }
  }
  return;
}

/*----------Main routine--------*/
void FSM(int *size, double *G, double *itt, double *velmod, double *out)
{
// Main variables
int Gnx, Gnz, Gnt, Gfull; // Sizes
int Gnx2,Gnz2; // Sizes +=2
double Gdx,Gdz,h; // Steps
double Got,Gdt; // Time info for rays
// Additional variables and arrays
int i,j; // Loop variables
double vmin, umax=1.e20; // Limit values for velocity and time
double **vel, **u; // Temporary 2D arrays for velocity and time

if(info!= 0) mexPrintf("Computing traveltimes with FSM...\n");
//Model sizes
Gnx = size[0]; Gnz = size[1];
Gfull = Gnx*Gnz; Gnx2=Gnx+2; Gnz2=Gnz+2; //+2 due to borders
//Steps
Gdx = G[6]; Gdz = G[8]; h = Gdx;
//Ray info
Got = G[9]; Gnt = (int)G[10]; Gdt = G[11];
vmin = minvalue(velmod, Gfull);
//Info
if(info!= 0) {
mexPrintf("Minimum velocity: Vmin=%f\n",vmin);
mexPrintf("Grid: Gnx=%d,Gnz=%d\n",Gnx,Gnz);
}
//Allocate temporary arrays
vel = new double*[Gnx2];
u = new double*[Gnx2];
for (i = 0; i < Gnx2; i++) {
  vel[i] = new double[Gnz2];
  u[i] = new double[Gnz2];
}
//Definition of velocity and traveltime
for (i=0;i<Gnx2;i++) {
  for (j=0;j<Gnz2;j++) {
    vel[i][j] = vmin;
      u[i][j] = umax;
  }
}
//Take velocity from input
for (i=0;i<Gnx;i++) {
  for (j=0;j<Gnz;j++) {
    vel[i+1][j+1] = velmod[ i + j*Gnx]; //shifting index by 1 to accomodate for borders
      u[i+1][j+1] =    itt[ i + j*Gnx];
  }
}

//Update times with FSM algorithm
fsmalgorithm(size,h,vel,u);

//Fill output array
for (i=0;i<Gnx;i++) {
  for (j=0;j<Gnz;j++) {
    out[ i + j*Gnx] = u[i+1][j+1]; //again shift by 1
  }
}
//Deallocate temporary arrays (to avoid memory leak)
for (i = 0; i < Gnx2; i++) {
  delete[]vel[i];
  delete[]u[i];
}
delete[] vel;
delete[] u;
//
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
int size[2];
double *G, *itt, *velmod, *tti;
//check input
if (nrhs!=3) mexErrMsgTxt("Three input arguments are required.");
if (nlhs!=1) mexErrMsgTxt("One output argument is required.");
//associate input pointers
G = (double *)mxGetPr(prhs[0]);
itt = (double *)mxGetPr(prhs[1]);
velmod = (double *)mxGetPr(prhs[2]);
//check if dimensions of velmod are correct
ndim = mxGetNumberOfDimensions(prhs[2]);
dims = mxGetDimensions(prhs[2]);
if (ndim!=2) mexErrMsgTxt("Number of dimensions of velmod array is not correct.");
if (info!= 0) mexPrintf("dims: dims[0]=%d, dims[1]=%d\n",dims[0],dims[1]);
size[0]=dims[0]; size[1]=dims[1];
// check if size of "itt" coinsides with size of "velmod" 
ndim = mxGetNumberOfDimensions(prhs[1]);
dims = mxGetDimensions(prhs[1]);
if (ndim!=2) mexErrMsgTxt("Number of dimensions of initial traveltime array (itt) is not correct.");
if (size[0]!=dims[0] || size[1]!=dims[1]){
    mexErrMsgTxt("Size of itt and velmod arrays should coincide");
  }
//create array for output
plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
//associate output pointer
tti = (double *)mxGetPr(plhs[0]);
// Call the master routine
if(info!= 0) mexPrintf("Calling FSM routine from Mex gateway...\n");
FSM(size, G, itt, velmod, tti);
if(info!= 0) mexPrintf("FSM finished.\n");

return;
}
