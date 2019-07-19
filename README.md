# MLIB
MATLAB has become a common tool in many research areas. MATLAB is a very user friendly programming language and suggests rich opportunities for development and research. For scientific computing, MATLAB includes a rich collection of built-in functions to solve a wide spectrum of mathematical problems. It's potential application can be further expanded by using optional toolboxes. As a result, a program written in MATLAB is relatively short and simple. Moreover, MATLAB offers powerful graphics capabilities, thus allowing easy to use visialization solutions to present the data.

Being an interpreted language, MATLAB executes instructions directly, without previously compiling a program into machine-language instructions. Hence, MATLAB is supposedly slow at loops compared to compiled languages like C, C++ or Fortran. However, functions written in the C/C++ or Fortran programming language, can be linked to MATLAB by using MEX library. MEX library allows to create MEX files from C/C++ or Fortran functions, that can be used and called directly from MATLAB as a built-in functions. Thus, a convenience of MATLAB can be combined with the high performance of C/C++ or Fortran languages. 

Reproducibility of research results is one of the key issues in the scientific community. MATLAB codes, together with supplementary content (comments, equations) and output can be published to the numerous formats, for example as PDF files. Such published documents can be used for sharing results with colleagues, teaching or demonstration as well as for generating comprehensible external documentation, thus making the research more reproducible. 

During the last years, MATLAB was my everyday tool for research, designing programs and testing ideas. I applied it for a variety of purposes, including:

Time imaging

Time lapse seismic

Traveltime tomography

Eikonal solvers, ray tracing

Localization of microseismic events

FD modelling, FWI

Static corrections

Virtual source method

I/O of X-Well, OBN, Microseismic, 2D/3D marine and field datasets. 

I found it helpful to collect the commonly used functions in the library, called MLIB. I also gathered Matlab-related open-source codes, that could be useful for seismic problems. In this script, I will demonstrate some tools of the MLIB library and several rarely used, but extremely useful features of MATLAB. In particular, I will demonstrate how to: 
add MLIB library read, write and analize SEGY/SU/SEG2 files use Eikonal solvers and ray tracing run 2D FD modeling perform quasi-Newton optimization method
I hope that MLIB library might be a usefull tool for fresh geophyscs students starting with MATLAB.
