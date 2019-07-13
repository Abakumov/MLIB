#!/bin/sh
#
# Abakumov Ivan
# 30.04.2015 Hamburg University
# e-mail: abakumov_ivan@mail.ru

# Take parameters corresponding to inline 390 from parameter cubes provided by Khawar
# Take only traces with CDP's: 

# if you want to make inline 390 only: 
mincdp=337659
maxcdp=338488

# if you want to make auxiliary volume around inline 390: 
#mincdp=325000
#maxcdp=351000

suwind < angle_sch.azimuth.su key=cdp min=$mincdp max=$maxcdp | susort > angle_sch_in390.azimuth.su cdp

suwind < angle_sch.dip.su key=cdp min=$mincdp max=$maxcdp | susort > angle_sch_in390.dip.su cdp

suwind < cmpini_sch.nip00.su key=cdp min=$mincdp max=$maxcdp | susort > cmpini_sch_in390.nip00.su cdp

suwind < cmpini_sch.nip10.su key=cdp min=$mincdp max=$maxcdp | susort > cmpini_sch_in390.nip10.su cdp

suwind < cmpini_sch.nip11.su key=cdp min=$mincdp max=$maxcdp | susort > cmpini_sch_in390.nip11.su cdp

suwind < crsstack_sch.stack.su key=cdp min=$mincdp max=$maxcdp | susort > crsstack_sch_in390.stack.su cdp

suwind < kn_sch.n00.su key=cdp min=$mincdp max=$maxcdp | susort > kn_sch_in390.n00.su cdp

suwind < kn_sch.n10.su key=cdp min=$mincdp max=$maxcdp | susort > kn_sch_in390.n10.su cdp

suwind < kn_sch.n11.su key=cdp min=$mincdp max=$maxcdp | susort > kn_sch_in390.n11.su cdp


#Finally copy manually results to:  
#/home/zmaw/u250128/Desktop/icrs_express/Schneeberg/product/khawar_old

