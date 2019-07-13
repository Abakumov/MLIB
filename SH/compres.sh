#!/bin/sh
#
# Abakumov Ivan
# 12.05.2015 Hamburg University
# e-mail: abakumov_ivan@mail.ru

# Compare results of crs and i-crs 

# crs (old)
file_old='/home/zmaw/u250128/Desktop/icrs_express/Schneeberg/product/test/test13/crsstack.cdp.stack.su'

# i-crs (new)
file_new='/home/zmaw/u250128/Desktop/icrs_express/Schneeberg/product/test/test14/icrsstack.cdp.stack.su'
#file_new='/scratch/local1/ivan/CRS_data/Schneeberg_3D_CRS_attributes/Khawar_old/crsstack_sch_x650.stack.su'
# crs (old)
#file_old='/home/zmaw/u250128/Desktop/icrs_express/Schneeberg/product/crs/test1/crsstack.cdp.coher.su'
# i-crs (new)
#file_new='/home/zmaw/u250128/Desktop/icrs_express/Schneeberg/product/crs/test7/icrsstack.cdp.coher.su'

# crs (old)
#file_old='/home/zmaw/u250128/Desktop/icrs_express/Schneeberg/product/crs/test1/crsstack.cdp.numt.su'
# i-crs (new)
#file_new='/home/zmaw/u250128/Desktop/icrs_express/Schneeberg/product/crs/test7/icrsstack.cdp.numt.su'




# STACK
suximage < $file_old wclip=-0.005 bclip=0.005 title=$file_old legend=1 cmap=rgb1 &
suximage < $file_new wclip=-0.005 bclip=0.005 title=$file_new legend=1 cmap=rgb1 &
sudiff $file_old $file_new | suximage wclip=-0.005 bclip=0.005 title=$file_new legend=1 cmap=rgb1 &

# COHER
#suximage < $file_old wclip=-0.000 bclip=0.02 title=$file_old legend=1 cmap=rgb1 &
#suximage < $file_new wclip=-0.000 bclip=0.02 title=$file_new legend=1 cmap=rgb1 &
#sudiff $file_old $file_new | suximage wclip=-0.005 bclip=0.005 title=$file_new legend=1 cmap=rgb1 &

# NUMT
#suximage < $file_old wclip=-0.000 bclip=5000 title=$file_old legend=1 cmap=rgb1 &
#suximage < $file_new wclip=-0.000 bclip=5000 title=$file_new legend=1 cmap=rgb1 &
#sudiff $file_old $file_new | suximage wclip=-10 bclip=10 title=$file_new legend=1 cmap=rgb1 &


