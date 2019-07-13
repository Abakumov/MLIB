#!/bin/sh
#
# Abakumov Ivan
# 15.05.2015 Hamburg University
# e-mail: abakumov_ivan@mail.ru

# Take parameters corresponding to xline 650 from parameter cubes provided by Khawar
# Take only traces with CDP's: 

# if you want to make xline 650 only: 
# cdp = 868*(inline-1) + xline; 

# Working folder is 
folder='/scratch/local1/ivan/CRS_data/Schneeberg_3D_CRS_attributes/Khawar_old'

# For each of parameter files (8 attributes + 1 stack)
afile=10;
while [ "$afile" -lt "11" ]
do 
  if [ "$afile" = "1" ]
    then
      filein='angle_sch.azimuth.su' 
      fileout='angle_sch_x650.azimuth.su'
  elif [ "$afile" = "2" ]
    then
      filein='angle_sch.dip.su'
      fileout='angle_sch_x650.dip.su'
  elif [ "$afile" = "3" ]
    then
      filein='cmpini_sch.nip00.su'
      fileout='cmpini_sch_x650.nip00.su'
  elif [ "$afile" = "4" ]
    then
      filein='cmpini_sch.nip10.su'
      fileout='cmpini_sch_x650.nip10.su'
  elif [ "$afile" = "5" ]
    then
      filein='cmpini_sch.nip11.su'
      fileout='cmpini_sch_x650.nip11.su'
  elif [ "$afile" = "6" ]
    then
      filein='kn_sch.n00.su'
      fileout='kn_sch_x650.n00.su'
  elif [ "$afile" = "7" ]
    then
      filein='kn_sch.n10.su'
      fileout='kn_sch_x650.n10.su'
  elif [ "$afile" = "8" ]
    then
      filein='kn_sch.n11.su'
      fileout='kn_sch_x650.n11.su'
  elif [ "$afile" = "9" ]
    then
      filein='crsstack_sch.stack.su'
      fileout='crsstack_sch_x650.stack.su'
  elif [ "$afile" = "10" ]
    then
      filein='crsstack_sch.numt.su'
      fileout='crsstack_sch_x650.numt.su'
  else 
    echo "Some strange error" 
  fi  

  fin=$folder/$filein
  fout=$folder/$fileout

  echo $fin
  echo $fout 


  # Cut xline 650
  
  inline=1
  while [ "$inline" -lt "713" ]
  do
    mincdp=`bc -l <<-END
       ($inline-1)*868+650	
    END`
 
    maxcdp=`bc -l <<-END
       ($inline-1)*868+650
    END`

    suwind < $fin key=cdp min=$mincdp max=$maxcdp > tmp$inline.su

    echo "Level $afile of 9: $inline out of 712  is done"

    inline=`expr $inline + 1`

  done

  #Part 10: Combine all files 
  cat tmp*.su > tmp
  susort <tmp cdp offset > $fout
  rm tmp*.su 
  rm tmp
  afile=`expr $afile + 1`
done




