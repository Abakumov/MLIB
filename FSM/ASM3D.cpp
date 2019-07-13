// FSM-type solver of adjoint eikonal equation
// version: 1.0
// Abakumov Ivan
// St. Petersburg University
// e-mail: abakumov_ivan@mail.ru
// 9th of July 2013
// Houston 
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
#define np 20 // Number of points around the shot where vel is "constant". default value 2 or 3
#define nploc 50
#define tmax 1e20

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////        AUXILIARY ROUTINES         /////////////////////////
///////////////////////////////////////////////////////////////////////////

int intround(double num) {
  return (num > 0.0) ? (int)floor(num + 0.5) : (int)ceil(num - 0.5);
}

///////////////////////////////////////////////////////////////////////////

double round(double num) {
  return (num > 0.0) ? floor(num + 0.5) : ceil(num - 0.5);
}

///////////////////////////////////////////////////////////////////////////

double minvalue(double* array, int size)
{
  static double val = array[0];
  for (int i = 1; i < size; i++) {
    val = val <= array[i] ? val : array[i];
  }
  return val;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//////  A D J O I N T   S T A T E   M E T H O D   /////////////////////////
///////////////////////////////////////////////////////////////////////////

double localasm(int i, int j, int k, double ***tti, double ***lambda, double ***DCR, double *G)
{
  double ap, am, app, apm, amp, amm;
  double bp, bm, bpp, bpm, bmp, bmm;
  double cp, cm, cpp, cpm, cmp, cmm;
  double coef, nlambda;
  double Gdx,Gdy,Gdz;         // Initial point and steps
  
  Gdx=G[6]; Gdy=G[7]; Gdz=G[8]; 

  ap = -(tti[i+1][j][k]-tti[i][j][k])/Gdx; app=(ap+fabs(ap))/2.; amp=(ap-fabs(ap))/2.;
  am = -(tti[i][j][k]-tti[i-1][j][k])/Gdx; apm=(am+fabs(am))/2.; amm=(am-fabs(am))/2.;
  bp = -(tti[i][j+1][k]-tti[i][j][k])/Gdy; bpp=(bp+fabs(bp))/2.; bmp=(bp-fabs(bp))/2.;
  bm = -(tti[i][j][k]-tti[i][j-1][k])/Gdy; bpm=(bm+fabs(bm))/2.; bmm=(bm-fabs(bm))/2.;
  cp = -(tti[i][j][k+1]-tti[i][j][k])/Gdz; cpp=(cp+fabs(cp))/2.; cmp=(cp-fabs(cp))/2.;
  cm = -(tti[i][j][k]-tti[i][j][k-1])/Gdz; cpm=(cm+fabs(cm))/2.; cmm=(cm-fabs(cm))/2.;
  
  coef = (app-amm)/Gdx + (bpp-bmm)/Gdy + (cpp-cmm)/Gdz;
  if (coef!=0){
  nlambda = (apm*lambda[i-1][j][k] - amp*lambda[i+1][j][k])/Gdx/coef + (bpm*lambda[i][j-1][k] - bmp*lambda[i][j+1][k])/Gdy/coef + (cpm*lambda[i][j][k-1] - cmp*lambda[k][j][k+1])/Gdz/coef + DCR[i][j][k]/coef; 
  }
  
  return(nlambda);
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
 
void Getlambda(double *G, double ***tti, double ***DCR, double ***lambda, int *MXYZ)
{

  int i,j,k,q;
  int MINX, MINY, MINZ, MAXX, MAXY, MAXZ;
  double Gdx,Gdy,Gdz;         // Initial point and steps
  
  MINX=MXYZ[0];
  MAXX=MXYZ[1];
  MINY=MXYZ[2];
  MAXY=MXYZ[3];
  MINZ=MXYZ[4];
  MAXZ=MXYZ[5];

  Gdx=G[6]; Gdy=G[7]; Gdz=G[8];  

  //Define initial values of lambda
  for (i=MINX;i<MAXX;i++){
    for (j=MINY;j<MAXY;j++){
      for (k=MINZ;k<MAXZ;k++){
        lambda[i][j][k]=DCR[i][j][k]/(-(tti[i][j][k+1]-tti[i][j][k])/Gdz);   //normal_up = {0, 0, -1}    (px*n1+py*n2+pz*n3)  up
      }
    }
  }
 
  //Iteration procedure for lambda
  for (q=0;q<2;q++){
      //#1
    for (i=MINX;i<MAXX;i++){
      for (j=MINY;j<MAXY;j++){
        for (k=MINZ;k<MAXZ;k++){
          lambda[i][j][k]=localasm(i,j,k,tti,lambda,DCR,G);
        }
      }                           
    }
      //#2
    for (i=MINX;i<MAXX;i++){
      for (j=MINY;j<MAXY;j++){
        for (k=MAXZ-1;k>MINZ-1;k--){
          lambda[i][j][k]=localasm(i,j,k,tti,lambda,DCR,G);
        }
      }                           
    }
      //#3
    for (i=MINX;i<MAXX;i++){
      for (j=MAXY-1;j>MINY-1;j--){
        for (k=MINZ;k<MAXZ;k++){
          lambda[i][j][k]=localasm(i,j,k,tti,lambda,DCR,G);
        }
      }                           
    }
      //#4
    for (i=MINX;i<MAXX;i++){
      for (j=MAXY-1;j>MINY-1;j--){
        for (k=MAXZ-1;k>MINZ-1;k--){
          lambda[i][j][k]=localasm(i,j,k,tti,lambda,DCR,G);
        }
      }                           
    }
      //#5
    for (i=MAXX-1;i>MINX-1;i--){
      for (j=MINY;j<MAXY;j++){
        for (k=MINZ;k<MAXZ;k++){
          lambda[i][j][k]=localasm(i,j,k,tti,lambda,DCR,G);
        }
      }                           
    }
      //#6
    for (i=MAXX-1;i>MINX-1;i--){
      for (j=MINY;j<MAXY;j++){
        for (k=MAXZ-1;k>MINZ-1;k--){
          lambda[i][j][k]=localasm(i,j,k,tti,lambda,DCR,G);
        }
      }                           
    }
      //#7
    for (i=MAXX-1;i>MINX-1;i--){
      for (j=MAXY-1;j>MINY-1;j--){
        for (k=MINZ;k<MAXZ;k++){
          lambda[i][j][k]=localasm(i,j,k,tti,lambda,DCR,G);
        }
      }                           
    }
      //#8
    for (i=MAXX-1;i>MINX-1;i--){
      for (j=MAXY-1;j>MINY-1;j--){
        for (k=MAXZ-1;k>MINZ-1;k--){
          lambda[i][j][k]=localasm(i,j,k,tti,lambda,DCR,G);
        }
      }                           
    }  
  }
  return;
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

/*----------Main routine--------*/
void ASM(double *G, double *traveltime, double *rightside, double *out)
{
    // Main variables
  int    MXYZ[6];
  int    Gnx, Gny, Gnz;       
  
    // Additional variables and arrays
  int i,j,k;                         // Loop variables
  double ***tti, ***dcr, ***lambda;  // Temporary 3D arrays 
  
  Gnx = (int)G[3]; 
  Gny = (int)G[4]; 
  Gnz = (int)G[5]; 
  
  MXYZ[0] = 1;
  MXYZ[1] = Gnx-1;
  MXYZ[2] = 1;
  MXYZ[3] = Gny-1;
  MXYZ[4] = 1;
  MXYZ[5] = Gnz-1;
  
    //Allocate temporary arrays
  tti = new double**[Gnx];
  dcr = new double**[Gnx];
  lambda = new double**[Gnx];
  for (i=0;i<Gnx;i++){
    tti[i] = new double*[Gny];
    dcr[i] = new double*[Gny];
    lambda[i] = new double*[Gny];
  }
  for (i=0;i<Gnx;i++){
    for (j=0;j<Gny;j++){
      tti[i][j] = new double[Gnz];
      dcr[i][j] = new double[Gnz];
      lambda[i][j] = new double[Gnz];
    }
  }
  for (i=0;i<Gnx;i++){
    for (j=0;j<Gny;j++){
      for (k=0;k<Gnz;k++){
        tti[i][j][k] = traveltime[i+j*Gnx+k*Gnx*Gny];  
        dcr[i][j][k] = rightside[i+j*Gnx+k*Gnx*Gny];  
        lambda[i][j][k] = 0;
      }
    }
  }
 
    // ASM algorithm
  Getlambda(G, tti, dcr, lambda, MXYZ);
  
    //Fill output array
  for (i=0;i<Gnx;i++){
    for (j=0;j<Gny;j++){
      for (k=0;k<Gnz;k++){
        out[i+j*Gnx+k*Gnx*Gny] = lambda[i][j][k]; 
      }
    }
  }
  
    //Deallocate temporary arrays (to avoid memory leak)  
  for (i=0;i<Gnx;i++){
    for (j=0;j<Gny;j++){
      delete[] tti[i][j];
      delete[] dcr[i][j];
      delete[] lambda[i][j];
    }
  }
  for (i=0;i<Gnx;i++){
    delete[] tti[i];
    delete[] dcr[i];
    delete[] lambda[i];
  }
  delete[] tti;
  delete[] dcr;
  delete[] lambda;

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*---- Gateway routine-----*/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    //declare variables
    const mwSize *dims, *rdims;
    mwSize ndim, rndim;
    double *G, *dcr, *tti;    
    double *lambda;
    //check input
    if (nrhs!=3) mexErrMsgTxt("Three input arguments are required.");  
    if (nlhs!=1) mexErrMsgTxt("One output argument is required.");    
    //associate input pointers
    G      = (double *)mxGetPr(prhs[0]);
    dcr = (double *)mxGetPr(prhs[1]);    
    tti = (double *)mxGetPr(prhs[2]);
    //check if dimensions of velmod are correct
    
    rndim = mxGetNumberOfDimensions(prhs[1]);
    rdims = mxGetDimensions(prhs[1]);   
    if (rndim!=3) mexErrMsgTxt("Number of dimensions of DCR array is not correct.");
    if (info!=0) mexPrintf("rdims: rdims[0]=%d, rdims[1]=%d, rdims[2]=%d\n",rdims[0],rdims[1],rdims[2]);    
    if (rdims[0]!=(int)G[3] || rdims[1]!=(int)G[4] || rdims[2]!=(int)G[5]) mexErrMsgTxt("Size of dcr array is not correctly defined in G-array");
    
    
    ndim = mxGetNumberOfDimensions(prhs[2]);
    dims = mxGetDimensions(prhs[2]);   
    if (ndim!=3) mexErrMsgTxt("Number of dimensions of tti array is not correct.");
    if (info!=0) mexPrintf("dims: dims[0]=%d, dims[1]=%d, dims[2]=%d\n",dims[0],dims[1],dims[2]);    
    if (dims[0]!=(int)G[3] || dims[1]!=(int)G[4] || dims[2]!=(int)G[5]) mexErrMsgTxt("Size of tti array is not correctly defined in G-array");
    //create array for output
    plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
    //associate output pointer
    lambda = (double *)mxGetPr(plhs[0]);    
    // Call the master routine
    if(info!= 0) mexPrintf("Calling ASM routine from Mex gateway...\n");
    ASM(G,tti,dcr,lambda);
    if(info!= 0) mexPrintf("ASM finished.\n");
    return;
}












