%% el2d make
% Compiles el2d.cpp code
clear all; close all; clc;


%% Compilation
tic
fprintf(1,'Compiling routine <el2d.cpp>...');

% LD_LIBRARY_PATH - Linux variable withpath to shared libraries 
getenv('LD_LIBRARY_PATH')

if(isunix)
    % Unix, gcc
   mex -v -largeArrayDims el2d.cpp COMPFLAGS="\$CFLAGS -fopenmp -Dinfo=0 -O3"...
       CXXFLAGS="\$CXXFLAGS -fopenmp -O3 -I/soft/FFTW/fftw-3.3.6-pl2/include"...
       LDFLAGS="\$LDFLAGS -fopenmp -L/soft/FFTW/fftw-3.3.6-pl2/lib -lfftw3f -O3";                % path to fftw library
else
    % Windows, MS Visual C++
    cmd_opt = ['-largeArrayDims '                        ...
        'el2d.cpp '                               ...
        'COMPFLAGS="$COMPFLAGS /openmp /O2 /I." ' ...
        'LINKFLAGS="$LINKFLAGS libfftw3f-3.lib"'];
    fid = fopen('mex_cmd_opt.rsp','w');
    fprintf(fid,'%s',cmd_opt);
    fclose(fid);
    mex @mex_cmd_opt.rsp;
end

fprintf(1,'Compiled. ');
toc

