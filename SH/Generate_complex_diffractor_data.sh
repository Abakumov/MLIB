#! /bin/sh

set -x

# Programm um komplexeres Diffraktormodell zu erstellen

susynlv nt=1001 dt=0.004 nxo=81 dxo=0.025 nxm=401 dxm=0.0125 v00=2.0 fpeak=30.0 \
ref="1:-1.0,3.0;6.0,3.0" ref="10:1.999,1.0;2.001,1.0" ref="10:2.999,2.0;3.001,2.0" ref="10:1.499,3.2;1.501,3.2" ref="10:1.999,2.0;2.001,2.0" ref="10:2.999,0.5;3.001,0.5" ref="10:3.499,1.5;3.501,1.5"\
smooth=1 > linse_klein.su
suaddnoise < linse_klein.su sn=5 > sn5.linse_klein.su

exit
dxo=0.025 fuer Offsets  dxo=0 fuer ZO
HOMOGEN
3 Diffraktoren  diffractor_3.su  sn5.diff_3.su
2 Diffraktoren  diffractor_2.su  sn5.diff_2.su
6 Diffraktoren 1 Reflektor mit und ohne Offset (_zo) sowie normal 20m langer 'Strich' und 'klein' 2m langer Strich

INHOMOGEN
3 Diffraktoren  inhomo_diffractor_3.su  sn5.inhomo_diff_3.su  gamma=0.5
3 Diffraktoren  lat_diffractor_3.su  sn5.lat_diff_3.su  dvdx=1.0

ref="10:1.99,1.0;2.01,1.0" ref="10:2.99,2.0;3.01,2.0" ref="10:1.49,3.2;1.51,3.2" ref="10:1.99,2.0;2.01,2.0" ref="10:2.99,0.5;3.01,0.5" ref="10:3.49,1.5;3.51,1.5" Linse 20m
