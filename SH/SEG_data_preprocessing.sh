#!/bin/csh
#
# Abakumov Ivan
# 17.03.2015 Hamburg University
# e-mail: abakumov_ivan@mail.ru

# Preprocessing of SEG C3 Salt WA dataset
# 1) automatic gain control
# 2) remove first arrivals
# 3) correct source to zero time    
# 4) hi cut filter 
# 5) add seismic noise 
# 
# Not applied 
# - threshhold amplitudes
# - statics for source delay
# 
# Set input and output filenames
#
set suinfile='SEG_C3WA_ffid_1-1200.su'
set suoutfile='SEG_C3WA_ffid_1-1200_proc.su' 
#
#set suinfile='SEG_C3WA_ffid_1201-2395.su'
#set suoutfile='SEG_C3WA_ffid_1201-2395_proc.su' 
#
########################################################
# Set parameters
#
# 1 Automatic gain control with first power of time  
set tpow=1
#
# 2 Remove first arrivals
# first arrivals are mixture of direct and reflected/refracted wave
# dt = offset/vel
# velocity in water is 1.5 km/sec
# velocity in segiments 1.7 km/sec
# effective velocity is chousen to be 1.65 km/sec
set vred=1.65
set tshift = 350 
set resampfac=1
# invert velocity km/s
set invred=`echo $vred $resampfac | awk '{print $2/$1}'`
#
# 3 correct source to zero time  
# addition data shift in ms (delay)  (88 msec in literature,  however I found that 124 is better) 
set add_tshift=124
#########################################################

# step 1, 2
sugain < $suinfile tpow=$tpow agc=1 wagc=0.5 |\
suchw key1=tstat key2=offset a=$tshift b=$invred | sustatic hdrs=1 sign=1  > tmp1.su

# step 2
suchw < tmp1.su key1=tstat key2=offset a=$tshift b=$invred | sustatic hdrs=1 sign=-1  > tmp2.su

# step 3, 4
suchw < tmp2.su key1=tstat a=$add_tshift b=0 | sustatic hdrs=1 sign=1 | sufilter f=1,5,20,30 > tmp3.su

# step 5 
suaddnoise < tmp3.su sn=10 > $suoutfile
#suximage < $suoutfile perc=99 title=$suoutfile&

exit
