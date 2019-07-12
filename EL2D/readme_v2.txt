README_v2: 28.04.2014
Updates:
since v1 (09.04.2014): 
-el2d.cpp has been improved with 4th order time derivative (fto global constant added)
-a bug has been corrected which produced a constant time shift by 3 time samples
-stype flag changed: 0 - explosion source, 1 - horizontal force; 2 - vertical force; 3 - shear source
-stability criterion has been updated accordingly
-----------------------------------------------------------------------------------------------------------------------
This is a short documentation on el2d MATLAB MEX/C++ function.
The function is based on original standalone C code created by Ekkehart Tessmer.
The function is designed to be compiled inside MATLAB using MEX either on UNIX(tested with gcc) or in WINDOWS(tested with Visual C++).
The function is intended for wavefield modelling in 2D elastic media with pseudospectral Fourier method.
For more information on the method see:
Kosloff, D., Reshef, M., & Loewenthal, D. (1984). 
Elastic wave calculations by the Fourier method. 
Bulletin of the Seismological Society of America, 74(3), 875-891.
http://www.bssaonline.org/content/74/3/875.short
-----------------------------------------------------------------------------------------------------------------------

Script example.m shows basic usage of the el2d function including compilation and plotting the result.
The el2d function is best suitable for shared-memory systems since 
it can be launched for several sources in parallel due to openmp multithreading.
-----------------------------------------------------------------------------------------------------------------------

The function accepts two modes of modelling:
- modelling with point source
- modelling with boundary condition as a source
First input argument of el2d function controls which mode is used, it could be 0 or 1. 
0 is modelling with point source, 
1 is modelling with boundary condition as a source
-----------------------------------------------------------------------------------------------------------------------

The fourth order time derivative scheme has been implemented (flag fto in el2d.cpp) and is used as default.
It has higher accuracy and allow bigger timesteps, but is slower than standard second order scheme.
In order to switch to second order scheme, change fto flag to 0.
-----------------------------------------------------------------------------------------------------------------------

Typical run of the el2d function in first modelling mode (shown in example):
[p, ux, uz] = el2d(0,rho,vp,vs,dx,dz,damp,ilx,ilz,stype,isx,isz,fpeak,dt,nt,wt);
Description of input arguments in that case:
rho: 2D array of model densities (in kg/m3)
vp: 2D array of model P-wave velocities (in m/s)
vs: 2D array of model S-wave velocities (in m/s)
dx: space sampling in horizontal direction (in m)
dz: space sampling in vertical direction (in m)
damp: 
if empty: no absorbing layer is applied at the edges of the model
if number: size of absorbing layer in grid nodes
if array of two numbers: first is size of absorbing layer, second is strength of damping (for 30 grid nodes optimal strength is 0.075)
if array with length more than 2: user defined damping coefficients, size of array correspond to the size of absorbing layer in grid nodes
ilx:
if 1D array: list of x grid nodes for data output
if 2D array: first row x grid nodes for data output, second row z grid nodes for data output (suitable for seismogram output)
ilz:
if ilx is 1D array: list of appropriate z grid nodes for data output (together ilx and ilz form a grid suitable for snapshot output)
otherwise: empty
stype: type of the source, 0: explosion source, 1: horizontal force, 2: vertical force, 3: shear source
isx: list of x grid node coordinates of sources
isz: list of z grid node coordinates of sources
fpeak: peak frequency of the source wavelet (required either for creation of source wavelet (if not defined by user) or for stability criteria estimation)
dt: time sampling (in s)
nt: number of time steps
wt: wavelet type
if number: 1: Gaussian 2: 1st derivative if Gaussian ,3: Ricker(2nd derivivative of Gaussian), 4: Seismic Ricker 
if array: user defined wavelet, must be of a size nt

Output can be from 1 to 7 output arrays, namely:
p: pressure
ux (optional): horizontal displacement
uz (optional): vertical displacement
uxx(optional): horizontal derivative of horizontal displacement
uzz(optional): vertical derivative of vertical displacement
uxz(optional): vertical derivative of horizontal displacement
uzx(optional): horizontal derivative of vertical displacement
if ilz is not empty, all arrays are of a size [nlx,nlz,nt,ns], where nlx is the size of ilx, nlz is the size of ilz 
if ilz is empty, all arrays are of a size [nlx,nt,ns], where nlx is the size of ilx, ns is the size of isx and isz
-----------------------------------------------------------------------------------------------------------------------

Typical run of the el2d function in second modelling mode (not shown in example):
[p, ux, uz] = el2d(1,rho,vp,vs,dx,dz,damp,ilx,ilz,btype,ibx,ibz,fpeak,dt,nt,bc1,bc2);
first 9 arguments has the same role
btype: type of boundary condition, 0: pressure bc, 1: horizontal displacement bc, 2: vertical displacement bc, 3: both vertical and horizontal displacement bc;
ibx: list of x grid node coordinates of points where boundary condition is injected
ibz: list of x grid node coordinates of points where boundary condition is injected
fpeak: proposed peak frequency of the signal (required for stability criteria estimation)
dt: time sampling (in s)
nt: number of time steps
bc1: boundary condition according to btype must be array of a size (nb,nt,ns), where nb is size of ibx and ibz and ns is number of shots for which the boundary condition was written
bc2: (is optional) vertical displacement boundary condition, should be present and of the same size as bc1 if only btype is 3.

Output is pretty much the same as in the first modelling mode.
-----------------------------------------------------------------------------------------------------------------------

The stability criteria are:
1. There should be a minimum of two grid points per shortest wavelength
2. Expression gamma/(dt*pi*vmax*sqrt(1./(dx^2)+1./(dz^2))) should be greater than 1, 
where gamma is equal to sqrt(12) for fourth order time derivative (default) and 2 for second order time derivative.
Both criteria are checked in example.m and also inside el2d.cpp.
If they are violated during the check inside el2d.cpp, then the warning will appear.
-----------------------------------------------------------------------------------------------------------------------

Absorbing boundary condition is implemented (with some modifications) according to
Cerjan, C., Kosloff, D., Kosloff, R., & Reshef, M. (1985). 
A nonreflecting boundary condition for discrete acoustic and elastic wave equations. 
Geophysics, 50(4), 705-708.
http://library.seg.org/doi/abs/10.1190/1.1441945
-----------------------------------------------------------------------------------------------------------------------

General requirements:
MATLAB (v2011 and later, but most likely earlier versions would work as well)
Special requirements under Linux:
gcc (>4.4)
FFTW library (fftw.org), usually it is in the list of available packages (for example, libfftw3-3 on Debian)
Special requirements under Windows:
Visual C++ 2010 or later
Special requirements under Windows 32bit:
FFTW libraries (fftw.org) for Windows 32bit (64bit libraries are provided with this readme)
-----------------------------------------------------------------------------------------------------------------------

Contents of the folder:
example.m: MATLAB example script
el2d.cpp: C++/MEX function itself
ricker.m: MATLAB function for ricker wavelet
makeColorMap: Auxiliary MATLAB function for colormap creation
Following files are required under Windows 64bit system only:
fftw3.h: FFTW library header file
libfftwf-3.dll: Dynamic-link FFTW single precision library for 64bit system
libfftwf-3.def: definition file for DLL
libfftwf-3.exp: export file for DLL
libfftwf-3.lib: lib file related to DLL library, to which the C++ is really linked
-----------------------------------------------------------------------------------------------------------------------

Author of this README:
Denis Anikiev (Saint Petersburg State University / University of Hamburg)
Questions, bug reports and feedback are welcomed: 
denis.anikiev@zmaw.de