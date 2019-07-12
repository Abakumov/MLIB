// Example of C++ code with OpenMP and MEX getaway routine
// Abakumov Ivan
// FU Berlin
// e-mail: abakumov_ivan@mail.ru
// 29th of January 2019
// Berlin

//c++  libraries
#include "mex.h"
#include <omp.h>
#include <iostream>
#include <cstdio>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <cstdlib>

//constants
#define info 0

using namespace std;
using std::cerr;
using std::endl;
using std::ofstream;

/* Shared Variables */
double minval[100]; 

int min_omp(double *x, int sizex, double *minx)
{
    for (int i=0; i<100; i++)
    {
        minval[i] = 1.e20;
    }
   
     #pragma omp parallel
    {
        
        int id = omp_get_thread_num();
        int total = omp_get_num_threads();
        int i, n, start, stop; 
        double mymin;        
        n = (int)sizex/total;           // size of x / number of threads 
        
        start = id*n;

        if ( id != (total-1) )
        {
           stop = start + n;
        }
        else
        {
           stop = sizex;
        }
        
        mymin = x[start];
        
        for (i = start+1; i < stop; i++ )
        {
            if ( x[i] < mymin )
                mymin = x[i];
        }
        
        minval[id] = mymin;
        
        printf("id %d total %d start %d stop %d mymin %f\n", id, total, start, stop, mymin);
    }

    minx[0] = minval[0];

    for (int i = 1; i < 100; i++)
    {
        if ( minval[i] < minx[0] )
            minx[0] = minval[i];
	}

    return 0;
}
  
   
///////////////////////////////////////////////////////////////////////////

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *x, *minx;
    int sizex;
    const mwSize *dims;
    mwSize ndim;
    // Set dimensions of output array
    mwSize ndimo = 2, dimo[2]; 
    dimo[0] = 1; dimo[1] = 1;
    
    //check input
    if (nrhs!=1) mexErrMsgTxt("One input argument is required.");
    if (nlhs!=1) mexErrMsgTxt("One output argument is required.");

    //associate input pointers
    x = (double *)mxGetPr(prhs[0]);
    
    // Find number of elements in x
    ndim = mxGetNumberOfDimensions(prhs[0]);
    dims = mxGetDimensions(prhs[0]);   
    if (ndim!=2) mexErrMsgTxt("Number of dimensions of 'x' array is not correct.");
    if (info!=0) mexPrintf("dims: dims[0]=%d, dims[1]=%d\n",dims[0],dims[1]);   
    sizex = int(dims[0]*dims[1]);
    
    //create array for output
    plhs[0] = mxCreateNumericArray(ndimo, dimo, mxDOUBLE_CLASS, mxREAL);
    
    //associate output pointer
    minx = (double *)mxGetPr(plhs[0]);

    min_omp(x,sizex,minx);
    
    return;
}