#!/bin/sh
#
# Abakumov Ivan
# 14.05.2015 Hamburg University
# e-mail: abakumov_ivan@mail.ru

# Make subset of Schneeberg data:
# Take only traces with CDP's: 
# X-Line 650, +-400 m,
# In-Line vary from 1 to 712 
# X-Line every 15 meters (400 m is approx 25 X-Lines)
# Required X-Line range: 625:675 
# Required In-Line range: 1:712

# Rule how the cdp number depend on In-Line and X-Line: 

# cdp = 868*(inline-1) + xline; 

# Input file:
filename='/scratch/local1/ivan/CRS_data/temp/outfile.su'

##############################################################################
# Part 1: Make two files out of one

iter=1
inline=1
dinline=384
while [ "$iter" -lt "3" ]
do

mincdp=`bc -l <<-END
   ($inline-1)*868+1	
END`
 
maxcdp=`bc -l <<-END
   ($inline-1+$dinline-1)*868+868
END`

suwind < $filename key=cdp min=$mincdp max=$maxcdp > IvanA$iter.su

echo "Level A $iter out of 2 is done"

inline=`expr $inline + $dinline`

iter=`expr $iter + 1`

done

##############################################################################
# Part 2: Make 4 files out of 2

iter=1
inline=1
dinline=192
while [ "$iter" -lt "5" ]
do

filenumber=`expr $iter + 1`
filenumber=`expr $filenumber / 2`


mincdp=`bc -l <<-END
   ($inline-1)*868+1	
END`
 
maxcdp=`bc -l <<-END
   ($inline-1+$dinline-1)*868+868
END`

suwind < IvanA$filenumber.su key=cdp min=$mincdp max=$maxcdp > IvanB$iter.su

echo "Level B $iter out of 4  is done"

inline=`expr $inline + $dinline`

iter=`expr $iter + 1`

done

rm IvanA*.su


##############################################################################
# Part 3: Make 8 files out of 4

iter=1
inline=1
dinline=96
while [ "$iter" -lt "9" ]
do

filenumber=`expr $iter + 1`
filenumber=`expr $filenumber / 2`

mincdp=`bc -l <<-END
   ($inline-1)*868+1	
END`
 
maxcdp=`bc -l <<-END
   ($inline-1+$dinline-1)*868+868
END`

suwind < IvanB$filenumber.su key=cdp min=$mincdp max=$maxcdp > IvanC$iter.su

echo "Level C $iter out of 8  is done"

inline=`expr $inline + $dinline`

iter=`expr $iter + 1`

done

rm IvanB*.su

##############################################################################
# Part 4: Make 16 files out of 8

iter=1
inline=1
dinline=48
while [ "$iter" -lt "16" ]
do

filenumber=`expr $iter + 1`
filenumber=`expr $filenumber / 2`

mincdp=`bc -l <<-END
   ($inline-1)*868+1	
END`
 
maxcdp=`bc -l <<-END
   ($inline-1+$dinline-1)*868+868
END`

suwind < IvanC$filenumber.su key=cdp min=$mincdp max=$maxcdp > IvanD$iter.su

echo "Level D $iter out of 15  is done"

inline=`expr $inline + $dinline`

iter=`expr $iter + 1`

done

rm IvanC*.su


##############################################################################
# Part 5: Make 30 files out of 15

iter=1
inline=1
dinline=24
while [ "$iter" -lt "31" ]
do

filenumber=`expr $iter + 1`
filenumber=`expr $filenumber / 2`

mincdp=`bc -l <<-END
   ($inline-1)*868+1	
END`
 
maxcdp=`bc -l <<-END
   ($inline-1+$dinline-1)*868+868
END`

suwind < IvanD$filenumber.su key=cdp min=$mincdp max=$maxcdp > IvanE$iter.su

echo "Level E $iter out of 30  is done"

inline=`expr $inline + $dinline`

iter=`expr $iter + 1`

done

rm IvanD*.su

##############################################################################
# Part 6: Make 60 files out of 30

iter=1
inline=1
dinline=12
while [ "$iter" -lt "61" ]
do

filenumber=`expr $iter + 1`
filenumber=`expr $filenumber / 2`

mincdp=`bc -l <<-END
   ($inline-1)*868+1	
END`
 
maxcdp=`bc -l <<-END
   ($inline-1+$dinline-1)*868+868
END`

suwind < IvanE$filenumber.su key=cdp min=$mincdp max=$maxcdp > IvanF$iter.su

echo "Level F $iter out of 60  is done"

inline=`expr $inline + $dinline`

iter=`expr $iter + 1`

done

rm IvanE*.su

##############################################################################
# Part 7: Make 120 files out of 60

iter=1
inline=1
dinline=6
while [ "$iter" -lt "120" ]
do

filenumber=`expr $iter + 1`
filenumber=`expr $filenumber / 2`

mincdp=`bc -l <<-END
   ($inline-1)*868+1	
END`
 
maxcdp=`bc -l <<-END
   ($inline-1+$dinline-1)*868+868
END`

suwind < IvanF$filenumber.su key=cdp min=$mincdp max=$maxcdp > IvanG$iter.su

echo "Level G $iter out of 119  is done"

inline=`expr $inline + $dinline`

iter=`expr $iter + 1`

done

rm IvanF*.su

##############################################################################
# Part 8: Make 238 files out of 119

iter=1
inline=1
dinline=3
while [ "$iter" -lt "239" ]
do

filenumber=`expr $iter + 1`
filenumber=`expr $filenumber / 2`

mincdp=`bc -l <<-END
   ($inline-1)*868+1	
END`
 
maxcdp=`bc -l <<-END
   ($inline-1+$dinline-1)*868+868
END`

suwind < IvanG$filenumber.su key=cdp min=$mincdp max=$maxcdp > IvanH$iter.su

echo "Level H $iter out of 239  is done"

inline=`expr $inline + $dinline`

iter=`expr $iter + 1`

done

rm IvanG*.su

##############################################################################
# Part 9: Make 712 files out of 238

iter=1
inline=1
dinline=1
while [ "$iter" -lt "713" ]
do

filenumber=`expr $iter + 2`
filenumber=`expr $filenumber / 3`

mincdp=`bc -l <<-END
   ($inline-1)*868+625	
END`
 
maxcdp=`bc -l <<-END
   ($inline-1+$dinline-1)*868+675
END`

suwind < IvanH$filenumber.su key=cdp min=$mincdp max=$maxcdp > tmp$iter.su

echo "Final level: $iter out of 712  is done"

inline=`expr $inline + $dinline`

iter=`expr $iter + 1`

done

rm IvanH*.su

##################################################################################
# Part 10: Combine all files 
cat tmp*.su > tmp
susort <tmp cdp offset > Schneeberg_xline_650_v2.su
rm tmp*.su 
rm tmp

exit
