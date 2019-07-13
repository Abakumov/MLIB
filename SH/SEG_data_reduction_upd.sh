#!/bin/sh
#
# Abakumov Ivan
# 30.03.2015 Hamburg University
# e-mail: abakumov_ivan@mail.ru

# Make subset of initial data:
# Take only traces with CDP's: 
# X = 8000, In-Line 200, +-500 m, 
# Y = 6000, X-Line 300, +-250 m, 
#X: 7500 - 8500  => In-Line 188-212
#Y: 5750 - 6250  => X-Line 288-312

#cdp = 1 + 188*900 + 1 = 169202
#cdp = 1 + 213*900 + 1 = 191702

# File SEG_C3WA_ffid_1-1200_proc.su contains In-Lines 93 - 196
# File SEG_C3WA_ffid_1201-2395_proc.su contains In-Lines 189 - 285

##############################################################################
# Part 1:  In-Lines 188 - 212
#cdp = 1 + 188*900 + 1 = 169202
#cdp = 1 + 213*900 + 1 = 191702

suwind < SEG_C3WA_ffid_1-1200_proc_upd.su key=cdp min=169202 max=191701 > tmp1.su
suwind < SEG_C3WA_ffid_1201-2395_proc_upd.su key=cdp min=169202 max=191701 > tmp2.su

###############################################################################
# Part 2: For In-Lines 93:187 take X-Lines 228-312

inline=93

while [ "$inline" -lt "188" ]
do

mincdp=`bc -l <<-END
   1+$inline*900+288	
END`
 
maxcdp=`bc -l <<-END
   1+$inline*900+312	
END`

suwind < SEG_C3WA_ffid_1-1200_proc_upd.su key=cdp min=$mincdp max=$maxcdp > tmp$inline.su

echo "Inline $inline is done"

inline=`expr $inline + 1`

done

###############################################################################
# Part 3: For In-Lines 213:285 take X-Lines 228-312

inline=213

while [ "$inline" -lt "286" ]
do

mincdp=`bc -l <<-END
   1+$inline*900+288	
END`
 
maxcdp=`bc -l <<-END
   1+$inline*900+312	
END`

suwind < SEG_C3WA_ffid_1201-2395_proc_upd.su key=cdp min=$mincdp max=$maxcdp > tmp$inline.su

echo "Inline $inline is done"

inline=`expr $inline + 1`

done

##################################################################################
# Part 4: Combine all files 
cat tmp*.su > tmp
susort <tmp cdp offset > SEG_C3WA_ffid_proc_red_upd.su

rm tmp*.su 










exit
