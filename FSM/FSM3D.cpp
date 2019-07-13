// Eikonal Equation solver by Fast Sweeping Method
// version: 1.6
// Abakumov Ivan
// St. Petersburg University
// e-mail: abakumov_ivan@mail.ru
// 23rd of October 2014
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
#define np 5 // Number of points around the shot where vel is "constant". default value 2 or 3
#define nploc 8
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


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
/////////////  F A S T   S W E E P I N G   M E T H O D   //////////////////
///////////////////////////////////////////////////////////////////////////

  //Exact traveltime routine for zone around source
double traveltime(double t0, double *x1, double *x2, double v)
{
  double t;
  t = t0 + sqrt(pow(x1[0]-x2[0],2)+pow(x1[1]-x2[1],2)+pow(x1[2]-x2[2],2))/v;
  return(t);
}

///////////////////////////////////////////////////////////////////////////

  // Routine for evaluation of time using FSM formulae
double xeval(double a, double b, double c, double d)
{
  double x,x1,x2,roots[3];
  roots[0] = a; roots[1] = b; roots[2] = c;
  int elements = sizeof(roots) / sizeof(roots[0]); 
  std::sort(roots, roots + elements); //sorting
  a = roots[0]; b = roots[1]; c = roots[2]; 
  x1 = a + d;
  x2 = (a + b + sqrt(2*pow(d,2) - pow((a-b),2)))/2;
  if  (x1 < b) {
    x = x1;
  } else if (x2 < c) {
    x = x2;
  } else {
    x = (a + b + c + sqrt(3*pow(d,2) + pow((a + b + c),2) - 3*(a*a + b*b + c*c)))/3;
  }
  return(x);
}

///////////////////////////////////////////////////////////////////////////

// FSM algorithm routine 
void fsmalgorithm(double *G, double ***vel, double ***u)
{
  double h = G[6];
  int i,j,k,q;                       // Loop variables
  int xmin,xmax,ymin,ymax,zmin,zmax; // Index limits          
  double a,b,c,d;                    // Temporary  
  xmin=1; xmax=(int)G[3]+1; 
  ymin=1; ymax=(int)G[4]+1; 
  zmin=1; zmax=(int)G[5]+1;
  
  for (q=0;q<2;q++) { //2 iterations of FSM 
    //1    
    for (i=xmin;i<xmax;i++) {
      for (j=ymin;j<ymax;j++) {
	    for (k=zmin;k<zmax;k++) {
	      a = min(u[i-1][j][k], u[i+1][j][k]);
	      b = min(u[i][j-1][k], u[i][j+1][k]);
	      c = min(u[i][j][k-1], u[i][j][k+1]);
	      d = h/vel[i][j][k];		    
	      u[i][j][k] = min( u[i][j][k], xeval(a,b,c,d) );
	    }
      }
    }
    //2
    for (i=xmin;i<xmax;i++) {
      for (j=ymin;j<ymax;j++) {
	    for (k=zmax-1;k>zmin-1;k--) {
	      a = min(u[i-1][j][k], u[i+1][j][k]);
	      b = min(u[i][j-1][k], u[i][j+1][k]);
	      c = min(u[i][j][k-1], u[i][j][k+1]);
	      d = h/vel[i][j][k];		    
	      u[i][j][k] = min( u[i][j][k], xeval(a,b,c,d) );
	    }
      }
    }
    //3    
    for (i=xmin;i<xmax;i++) {
      for (j=ymax-1;j>ymin-1;j--) {
	    for (k=zmin;k<zmax;k++) {
	      a = min(u[i-1][j][k], u[i+1][j][k]);
	      b = min(u[i][j-1][k], u[i][j+1][k]);
	      c = min(u[i][j][k-1], u[i][j][k+1]);
	      d = h/vel[i][j][k];		    
	      u[i][j][k] = min( u[i][j][k], xeval(a,b,c,d) );
	    }
      }
    }
    //4
    for (i=xmin;i<xmax;i++) {
      for (j=ymax-1;j>ymin-1;j--) {
	    for (k=zmax-1;k>zmin-1;k--) {
	      a = min(u[i-1][j][k], u[i+1][j][k]);
	      b = min(u[i][j-1][k], u[i][j+1][k]);
	      c = min(u[i][j][k-1], u[i][j][k+1]);
	      d = h/vel[i][j][k];		    
	      u[i][j][k] = min( u[i][j][k], xeval(a,b,c,d) );
	    }
      }
    }
    //5
    for (i=xmax-1;i>xmin-1;i--) {
      for (j=ymin;j<ymax;j++) {
	    for (k=zmin;k<zmax;k++) {
	      a = min(u[i-1][j][k], u[i+1][j][k]);
	      b = min(u[i][j-1][k], u[i][j+1][k]);
	      c = min(u[i][j][k-1], u[i][j][k+1]);
	      d = h/vel[i][j][k];		    
	      u[i][j][k] = min( u[i][j][k], xeval(a,b,c,d) );
	    }
      }
    }
    //6
    for (i=xmax-1;i>xmin-1;i--) {
      for (j=ymin;j<ymax;j++) {
	    for (k=zmax-1;k>zmin-1;k--) {
	      a = min(u[i-1][j][k], u[i+1][j][k]);
	      b = min(u[i][j-1][k], u[i][j+1][k]);
	      c = min(u[i][j][k-1], u[i][j][k+1]);
	      d = h/vel[i][j][k];		    
	      u[i][j][k] = min( u[i][j][k], xeval(a,b,c,d) );
	    }
      }
    }
    //7
    for (i=xmax-1;i>xmin-1;i--) {
      for (j=ymax-1;j>ymin-1;j--) {
	    for (k=zmin;k<zmax;k++) {
	      a = min(u[i-1][j][k], u[i+1][j][k]);
	      b = min(u[i][j-1][k], u[i][j+1][k]);
	      c = min(u[i][j][k-1], u[i][j][k+1]);
	      d = h/vel[i][j][k];		    
	      u[i][j][k] = min( u[i][j][k], xeval(a,b,c,d) );
	    }
      }
    }
    //8
    for (i=xmax-1;i>xmin-1;i--) {
      for (j=ymax-1;j>ymin-1;j--) {
	    for (k=zmax-1;k>zmin-1;k--) {
	      a = min(u[i-1][j][k], u[i+1][j][k]);
	      b = min(u[i][j-1][k], u[i][j+1][k]);
	      c = min(u[i][j][k-1], u[i][j][k+1]);
	      d = h/vel[i][j][k];		    
	      u[i][j][k] = min( u[i][j][k], xeval(a,b,c,d) );
	    }
      }
    }
  }
  return;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void localtimegrid(double *G, double *S, double ***vel, double ***tti)      
{
    int Gnx,Gny,Gnz;
    int Gsx,Gsy,Gsz;
    int i,j,k;
    int lx, ly, lz;
    int rx, ry, rz;
    int im, jm, km;
    int ip, jp, kp;
    int ii, jj, kk;
    int locx, locy, locz; 
    int locGnx, locGny, locGnz;
    int locGsx, locGsy, locGsz;
    int iloc, jloc, kloc;
    int locside;
  
        
    double fragm, fragm3;
    
    double locG[12];
    double Gdx, Gdy, Gdz;
    double Gox, Goy, Goz;
    double ***locvel;
    double ***loctti;
    double locGox, locGoy, locGoz;
    double locGdx, locGdy, locGdz;
    double veleff;
    double loc[3];
    double Got=0;
    
    Gox = G[0];      Goy = G[1];      Goz = G[2];
    Gnx =(int)G[3];  Gny = (int)G[4]; Gnz = (int)G[5];
    Gdx = G[6];      Gdy = G[7];      Gdz = G[8]; 
    
    if(info!= 0) {
        mexPrintf("Gnx=%d, Gny=%d, Gnz=%d\n",Gnx,Gny,Gnz);
    }
    
    if(S[0]<Gox-Gdx || S[0]>(Gox+Gdx*Gnx) || S[1]<Goy-Gdy || S[1]>(Goy+Gdy*Gny) || S[2]<Goz-Gdz || S[2]>(Goz+Gdz*Gnz)){
      mexErrMsgTxt("Shot coordinates out of range");
    }
    Gsx=intround((S[0]-Gox)/Gdx)+1;         // values on the new grid!
    Gsy=intround((S[1]-Goy)/Gdy)+1;
    Gsz=intround((S[2]-Goz)/Gdz)+1;
    
    if(info!= 0) {
        mexPrintf("Gsx=%d, Gsy=%d, Gsz=%d\n",Gsx,Gsy,Gsz);
    }  
    // example
    lx = max(Gsx-np, 0);        
    ly = max(Gsy-np, 0);        
    lz = max(Gsz-np, 0);        
    rx = min(Gsx+np, Gnx+1);    
    ry = min(Gsy+np, Gny+1);
    rz = min(Gsz+np, Gnz+1);    
                                
    fragm = 4.0;               
    fragm3 = fragm*fragm*fragm;
    locx = (rx-lx)*(int)fragm+1;         
    locy = (ry-ly)*(int)fragm+1;
    locz = (rz-lz)*(int)fragm+1;     
    
    if(info!= 0) {
        mexPrintf("lx=%d, rx=%d, ly=%d, ry=%d, lz=%d, rz=%d\n",lx,rx,ly,ry,lz,rz);
        mexPrintf("locx=%d, locy=%d, locz=%d\n",locx,locy,locz);
    }
    
    //Step 1.
    locvel = new double**[locx];
    loctti = new double**[locx];
    for (i=0;i<locx;i++){
      locvel[i] = new double*[locy];
      loctti[i] = new double*[locy];
    }           
    for (i=0;i<locx;i++){
      for (j=0;j<locy;j++){
        locvel[i][j] = new double[locz];
        loctti[i][j] = new double[locz];
      }
    }
    
    for (i=0;i<locx;i++){
      for (j=0;j<locy;j++){
        for (k=0;k<locz;k++){  
            // define loctti
          loctti[i][j][k] = tmax;
            // and interpolate velocity model 
                                       
          im = lx + (int)floor(i/fragm); 
          jm = ly + (int)floor(j/fragm); 
          km = lz + (int)floor(k/fragm); 
          ip = lx + (int)ceil(i/fragm);  
          jp = ly + (int)ceil(j/fragm);          
          kp = lz + (int)ceil(k/fragm);
          
          ii = i - (im-lx)*(int)fragm;
          jj = j - (jm-ly)*(int)fragm;
          kk = k - (km-lz)*(int)fragm;
        
          locvel[i][j][k] = vel[im][jm][km]/fragm3*(fragm-ii)*(fragm-jj)*(fragm-kk)+vel[im][jm][kp]/fragm3*(fragm-ii)*(fragm-jj)*kk+vel[im][jp][km]/fragm3*(fragm-ii)*jj*(fragm-kk)+vel[im][jp][kp]/fragm3*(fragm-ii)*jj*kk+vel[ip][jm][km]/fragm3*ii*(fragm-jj)*(fragm-kk)+vel[ip][jm][kp]/fragm3*ii*(fragm-jj)*kk+vel[ip][jp][km]/fragm3*ii*jj*(fragm-kk)+vel[ip][jp][kp]/fragm3*ii*jj*kk;
        }
      }
    }
    
   
    locGox = Gox + (lx-1)*Gdx;      // 0 (old) <-> Gox - Gdx (old)
    locGoy = Goy + (ly-1)*Gdy;      // 0 (old) <-> Gox - Gdx (old)
    locGoz = Goz + (lz-1)*Gdz;      // 1 (old) <-> Gox  (old)
                                    
    locGdx = Gdx/fragm;
    locGdy = Gdy/fragm;
    locGdz = Gdz/fragm;
    
    locGnx = locx;
    locGny = locy;
    locGnz = locz;
    
    locGsx=intround((S[0]-locGox)/locGdx)+1;
    locGsy=intround((S[1]-locGoy)/locGdy)+1;
    locGsz=intround((S[2]-locGoz)/locGdz)+1;
    
    locG[0] = locGox;
    locG[1] = locGoy;
    locG[2] = locGoz;
    locG[3] = (double)(locGnx-2);       // see how G is defined in FSM algorithm
    locG[4] = (double)(locGny-2);
    locG[5] = (double)(locGnz-2);
    locG[6] = locGdx;
    locG[7] = locGdy;
    locG[8] = locGdz;
    
    locside = 2*nploc+1;
    
    iloc=locGsx-nploc-1;
    jloc=locGsy-nploc-1;
    kloc=locGsz-nploc-1;
    
    for (i=max(0,iloc);i<min(locGnx,iloc+locside);i++) {
      for (j=max(0,jloc);j<min(locGny,jloc+locside);j++) {
        for (k=max(0,kloc);k<min(locGnz,kloc+locside);k++) {
          loc[0] = locGox + i*locGdx;  
          loc[1] = locGoy + j*locGdy;
          loc[2] = locGoz + k*locGdz; 
          veleff = 2*locvel[i][j][k]*locvel[locGsx][locGsy][locGsz]/(vel[i][j][k]+vel[locGsx][locGsy][locGsz]);
          loctti[i][j][k]=traveltime(Got,S,loc,veleff); 
        }
      }    
    }
      
    // start FSM algorithm    
    fsmalgorithm(locG, locvel, loctti);
    
    // assign values from local grid to normal grid
    for (i=0;i<locx;){
      im = (int)floor(i/fragm);   
      for (j=0;j<locy;){
        jm = (int)floor(j/fragm);  
        for (k=0;k<locz;){  
          km = (int)floor(k/fragm); 
          tti[im+lx][jm+ly][km+lz] = loctti[i][j][k];
          k = k +(int)fragm;
        }
        j = j + (int)fragm;
      }
      i = i + (int)fragm;
    }
       
    // delete locvel and loctti
     for (i=0;i<locx;i++){
       for (j=0;j<locy;j++){
         delete[]locvel[i][j];
         delete[]loctti[i][j];
       }
     }
     for (i=0;i<locx;i++){
       delete[]locvel[i];
       delete[]loctti[i];
     }
     delete[] locvel;
     delete[] loctti;
    
     return;
}


/*----------Main routine--------*/
void FSM(double *G, double *shot, double *velmod, double *out)
{
    // Main variables
  int Gnx, Gny, Gnz, Gnt, Gfull;     // Sizes
  int Gnx2, Gny2, Gnz2;              // Sizes +=2
  double Gdx, Gdy, Gdz, step;              // Steps
  
    // Additional variables and arrays
  int i,j,k;                         // Loop variables
  double vmin;                       // Limit values for velocity and time
  double ***vel, ***u;               // Temporary 3D arrays for velocity and time
  
  if(info!= 0) mexPrintf("Computing traveltimes with FSM...\n");
    
    //Model sizes  
  Gnx = (int)G[3]; Gny = (int)G[4]; Gnz = (int)G[5]; 
  Gfull = Gnx*Gny*Gnz; 
  Gnx2=Gnx+2; 
  Gny2=Gny+2; 
  Gnz2=Gnz+2; //+2 due to borders

    //Steps
  Gdx = G[6]; Gdy = G[7]; Gdz = G[8]; step = Gdx;

  vmin = minvalue(velmod, Gfull);
    
    //Info
  if(info!= 0) {
  mexPrintf("Minimum velocity: Vmin=%f\n",vmin);
  mexPrintf("Grid:             Gnx=%d,Gny=%d,Gnz=%d\n",Gnx,Gny,Gnz);
  mexPrintf("Source position:  sx=%f,sy=%f,sz=%f\n",shot[0],shot[1],shot[2]);
  }
  
    //Allocate temporary arrays
  vel = new double**[Gnx2];
  u   = new double**[Gnx2];
  for (i = 0; i < Gnx2; i++) {
    vel[i] = new double*[Gny2];
    u[i]   = new double*[Gny2];
  }
  for (i = 0; i < Gnx2; i++) {
    for (j = 0; j < Gny2; j++) {
      vel[i][j] = new double[Gnz2];
      u[i][j]   = new double[Gnz2];
    }
  }

    //Definition of velocity and traveltime  
  for (i=0;i<Gnx2;i++) {
    for (j=0;j<Gny2;j++) {
      for (k=0;k<Gnz2;k++) {
    	vel[i][j][k] = vmin;
	    u[i][j][k]   = tmax;
      }
    }
  }
  
    //Take velocity from input      
  for (i=0;i<Gnx;i++) {
    for (j=0;j<Gny;j++) {
      for (k=0;k<Gnz;k++) {
        vel[i+1][j+1][k+1] = velmod[i+j*Gnx+k*Gnx*Gny];  //shifting index by 1 to accomodate for borders        
      }
    }
  }
    //Initialize times in source zone manually
  localtimegrid(G,shot,vel,u);
  
    //Update times with FSM algorithm
  fsmalgorithm(G,vel,u);    
  
    //Fill output array
  for (i=0;i<Gnx;i++) {
    for (j=0;j<Gny;j++) {
      for (k=0;k<Gnz;k++) {
        out[i+j*Gnx+k*Gnx*Gny] = u[i+1][j+1][k+1];  //again shift by 1	           
      }
    }
  }
  
//   mexPrintf("Save Traveltimes\n");
//   FILE *fl1=NULL;
//   fl1 = fopen("tmpTTp.dat", "wb+");
//   fwrite(out,sizeof(double),Gfull,fl1);
//   fclose(fl1);   

    //Deallocate temporary arrays (to avoid memory leak)  
  for (i = 0; i < Gnx2; i++) {
    for (j = 0; j < Gny2; j++) {
      delete[]vel[i][j];
      delete[]u[i][j];
    }
  }
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
    double *G, *shot, *velmod, *tti;    
    //check input
    if (nrhs!=3) mexErrMsgTxt("Three input arguments are required.");  
    if (nlhs!=1) mexErrMsgTxt("One output argument is required.");    
    //associate input pointers
    G      = (double *)mxGetPr(prhs[0]);
    shot   = (double *)mxGetPr(prhs[1]);    
    velmod = (double *)mxGetPr(prhs[2]);
    //check if dimensions of velmod are correct
    ndim = mxGetNumberOfDimensions(prhs[2]);
    dims = mxGetDimensions(prhs[2]);   
    if (ndim!=3) mexErrMsgTxt("Number of dimensions of velmod array is not correct.");
    if (info!=0) mexPrintf("dims: dims[0]=%d, dims[1]=%d, dims[2]=%d\n",dims[0],dims[1],dims[2]);    
    if (dims[0]!=(int)G[3] || dims[1]!=(int)G[4] || dims[2]!=(int)G[5]) mexErrMsgTxt("Size of velocity array is not correctly defined in G-array");
    //create array for output
    plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
    //associate output pointer
    tti = (double *)mxGetPr(plhs[0]);    
    // Call the master routine
    if(info!= 0) mexPrintf("Calling FSM routine from Mex gateway...\n");
    FSM(G, shot, velmod, tti);
    if(info!= 0) mexPrintf("FSM finished.\n");
    return;
}
