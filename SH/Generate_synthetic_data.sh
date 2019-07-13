#!/bin/sh
#
# Abakumov Ivan
# 22.06.2015 Hamburg University
# e-mail: abakumov_ivan@mail.ru

# This script generates simple PP and PS seismograms for different types of reflectors and different velocity models

########################################################################################################
################## What kind of reflector do you like?

#ref='plane'				# plane reflector, 1km depth
#ref='r10'				# curved reflector, 10 m radius, top 1km, plane part 1.01 km depth 
#ref='r100'				# curved reflector, 100 m radius, top 1km, plane part 1.1 km depth
ref='r1000'				# curved reflector, 1km radius, top 1 km, plane part 2 km depth
#ref='r10000'				# curved reflector, 10 km radius, top 1km

# Note: it is possible to impruve reflectors (dots are not optimal and in some cases not simmetric). 

################## What kind of waves we will model? 

wave=PP					# PP wave
#wave=PS				# PS wave

################# Gradient needed?

dvdz=0.0				# no gradient
#dvdz=0.5
#dvdz=1.0
#dvdz=1.5 				# very strong gradient

################ Set other properties
# nt                 number of time samples
# dt                 time sampling interval (sec)
# nxo                number of source-receiver offsets
# dxo                offset sampling interval (kilounits)
# nxm                number of midpoints (see notes below)
# dxm                midpoint sampling interval (kilounits)
# x0=0.0                 distance x at which v00 is specified
# z0=0.0                 depth z at which v00 is specified
# v00=2.0                velocity at x0,z0 (kilounits/sec)
# dvdx=0.0               derivative of velocity with distance x (dv/dx)
# dvdz=0.0               derivative of velocity with depth z (dv/dz)
# fpeak=0.2/dt           peak frequency of symmetric Ricker wavelet (Hz)
# ref="1:1,2;4,2"        reflector(s):  "amplitude:x1,z1;x2,z2;x3,z3;..."
# smooth=0               =1 for smooth (piecewise cubic spline) reflectors
# er=0                   =1 for exploding reflector amplitudes
# ls=0                   =1 for line source; default is point source

prop='nt=1001 dt=0.004 nxo=81 dxo=0.025 nxm=201 dxm=0.025 fpeak=50.0 v00=2.0'
gamma=0.6667

########################################################################################################
################################### No changes below this line.... #####################################

if [ "$ref" = "plane" ]
    then
	reflector='-1, 1; 6, 1'
elif [ "$ref" = "r10" ]
    then
	reflector='-1,1.010000; 2.490000,1.010000; 2.490500,1.006878; 2.491000,1.005641; 2.491500,1.004732; 2.492000,1.004000; 2.492500,1.003386; 2.493000,1.002859; 2.493500,1.002401; 2.494000,1.002000; 2.494500,1.001648; 2.495000,1.001340; 2.495500,1.001070; 2.496000,1.000835; 2.496500,1.000633; 2.497000,1.000461; 2.497500,1.000318; 2.498000,1.000202; 2.498500,1.000113; 2.499000,1.000050; 2.499500,1.000013; 2.500000,1.000000; 2.500500,1.000013; 2.501000,1.000050; 2.501500,1.000113; 2.502000,1.000202; 2.502500,1.000318; 2.503000,1.000461; 2.503500,1.000633; 2.504000,1.000835; 2.504500,1.001070; 2.505000,1.001340; 2.505500,1.001648; 2.506000,1.002000; 2.506500,1.002401; 2.507000,1.002859; 2.507500,1.003386; 2.508000,1.004000; 2.508500,1.004732; 2.509000,1.005641; 2.509500,1.006878; 2.510000,1.010000; 6.000000,1.010000'
elif [ "$ref" = "r100" ]
    then
	#reflector='-1,1.100000;2.400000,1.100000;2.425000,1.033856;2.450000,1.013397;2.475000,1.003175;2.500000,1.000000;2.525000,1.003175;2.550000,1.013397;2.575000,1.033856;2.600000,1.100000;6,1.100000;'
	reflector='2.400000,1.100000;2.425000,1.033856;2.450000,1.013397;2.475000,1.003175;2.500000,1.000000;2.525000,1.003175;2.550000,1.013397;2.575000,1.033856;2.600000,1.100000;'
elif [ "$ref" = "r1000" ]
    then
	reflector='-1, 2; 0.85, 2; 0.9, 2; 0.95, 2; 1.0, 2; 1.05, 2; 1.1, 2; 1.15, 2; 1.2, 2; 1.25, 2; 1.3, 2; 1.35, 2; 1.4, 2; 1.45, 2; 1.5, 2; 1.501, 1.96; 1.504, 1.91; 1.507, 1.88; 1.513, 1.84; 1.519, 1.80; 1.526, 1.77; 1.538, 1.73; 1.550, 1.68; 1.563, 1.65; 1.58, 1.62; 1.588, 1.59; 1.601, 1.56; 1.626, 1.51; 1.651, 1.47; 1.676, 1.43; 1.701, 1.40; 1.726, 1.37; 1.751, 1.34; 1.776, 1.31; 1.826, 1.26; 1.876, 1.22; 1.926, 1.18; 1.976, 1.15; 2.0258, 1.12; 2.075, 1.09; 2.126, 1.07; 2.176, 1.05; 2.226, 1.04; 2.276, 1.025; 2.326,  1.02;  2.376,  1.01;  2.426, 1.003; 2.5, 1.00; 2.576,  1.003;  2.626,  1.01; 2.676, 1.02; 2.726,  1.03; 2.776, 1.04; 2.826, 1.05; 2.876, 1.07; 2.926, 1.09; 2.976, 1.12; 3.026, 1.15; 3.076, 1.18; 3.126, 1.22; 3.176, 1.26; 3.226, 1.31; 3.251, 1.34; 3.276, 1.37; 3.301, 1.40;  3.326, 1.44; 3.351, 1.47; 3.376, 1.52; 3.401, 1.57; 3.413, 1.59; 3.426, 1.62; 3.438, 1.65; 3.451, 1.69; 3.463, 1.73; 3.475, 1.781; 3.482, 1.81; 3.488, 1.85; 3.494, 1.90; 3.497, 1.93; 3.499, 1.96; 3.55, 2; 6, 2'
elif [ "$ref" = "r10000" ]
    then
	reflector='-1.000000,1.632503;0.000000,1.317542;1.000000,1.113140;2.000000,1.012508;3.000000,1.012508;4.000000,1.113140;5.000000,1.317542;6.000000,1.632503'
else 
    echo "Reflector not properly defined" 
fi 

if [ "$wave" = "PP" ]
    then
	code='susynlv'
elif [ "$wave" = "PS" ]
    then
	code='susynlvcw'
else 
    echo "Modeling code is not defined" 
fi 

filename=$wave.grad.$dvdz.$ref.su

echo $code is used
echo Output filename - $filename 

$code > $filename $prop dvdz=$dvdz gamma=$gamma ref="1: $reflector"
