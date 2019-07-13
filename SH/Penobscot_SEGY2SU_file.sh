#!/bin/sh
#
# Abakumov Ivan
# 18.06.2015 Hamburg University
# e-mail: abakumov_ivan@mail.ru

# Penobscot data
# original files in folder SEGY
#-rwxrwxrwx 1 u250115 u250 3811416128 Jun 18 14:11 3D_gathers_pstm_nmo_X1001.sgy
#-rwxrwxrwx 1 u250115 u250 6170492988 Jun 18 14:19 3D_gathers_pstm_nmo_X1041.sgy
#-rwxrwxrwx 1 u250115 u250 7309292440 Jun 18 14:30 3D_gathers_pstm_nmo_X1081.sgy
#-rwxrwxrwx 1 u250115 u250 7322935580 Jun 18 14:40 3D_gathers_pstm_nmo_X1121.sgy
#-rwxrwxrwx 1 u250115 u250 7322948068 Jun 18 14:51 3D_gathers_pstm_nmo_X1161.sgy
#-rwxrwxrwx 1 u250115 u250 7322873140 Jun 18 14:05 3D_gathers_pstm_nmo_X1201.sgy
#-rwxrwxrwx 1 u250115 u250 7321355848 Jun 18 15:01 3D_gathers_pstm_nmo_X1241.sgy
#-rwxrwxrwx 1 u250115 u250 7323584956 Jun 18 15:12 3D_gathers_pstm_nmo_X1281.sgy
#-rwxrwxrwx 1 u250115 u250 7323678616 Jun 18 15:22 3D_gathers_pstm_nmo_X1321.sgy
#-rwxrwxrwx 1 u250115 u250 7320032120 Jun 18 15:33 3D_gathers_pstm_nmo_X1361.sgy
#-rwxrwxrwx 1 u250115 u250 7317821744 Jun 18 15:43 3D_gathers_pstm_nmo_X1401.sgy
#-rwxrwxrwx 1 u250115 u250 7313544604 Jun 18 15:54 3D_gathers_pstm_nmo_X1441.sgy
#-rwxrwxrwx 1 u250115 u250 7302080620 Jun 18 16:04 3D_gathers_pstm_nmo_X1481.sgy
#-rwxrwxrwx 1 u250115 u250 6179209612 Jun 18 16:13 3D_gathers_pstm_nmo_X1521.sgy
#-rwxrwxrwx 1 u250115 u250 3762506876 Jun 18 16:18 3D_gathers_pstm_nmo_X1561.sgy

# Goal of this script - read files and create one big su file and one big file with headers
# Save results in folder SU

segypath='/scratch/local1/ivan/CRS_data/Penobscot_3D_data/3D_gathers_pstm_nmo_sgy/SEGY'
supath='/scratch/local1/ivan/CRS_data/Penobscot_3D_data/3D_gathers_pstm_nmo_sgy/SU'
afile=1;
while [ "$afile" -lt "16" ]
do 
  if [ "$afile" = "1" ]
    then
	file='3D_gathers_pstm_nmo_X1001'      
  elif [ "$afile" = "2" ]
    then
	file='3D_gathers_pstm_nmo_X1041'      
  elif [ "$afile" = "3" ]
    then
	file='3D_gathers_pstm_nmo_X1081'      
  elif [ "$afile" = "4" ]
    then
	file='3D_gathers_pstm_nmo_X1121'      
  elif [ "$afile" = "5" ]
    then
	file='3D_gathers_pstm_nmo_X1161'      
  elif [ "$afile" = "6" ]
    then
	file='3D_gathers_pstm_nmo_X1201'      
  elif [ "$afile" = "7" ]
    then
	file='3D_gathers_pstm_nmo_X1241'      
  elif [ "$afile" = "8" ]
    then
	file='3D_gathers_pstm_nmo_X1281'      
  elif [ "$afile" = "9" ]
    then
	file='3D_gathers_pstm_nmo_X1321'      
  elif [ "$afile" = "10" ]
    then
	file='3D_gathers_pstm_nmo_X1361'      
  elif [ "$afile" = "11" ]
    then
	file='3D_gathers_pstm_nmo_X1401'      
  elif [ "$afile" = "12" ]
    then
	file='3D_gathers_pstm_nmo_X1441'      
  elif [ "$afile" = "13" ]
    then
	file='3D_gathers_pstm_nmo_X1481'      
  elif [ "$afile" = "14" ]
    then
	file='3D_gathers_pstm_nmo_X1521'      
  elif [ "$afile" = "15" ]
    then
	file='3D_gathers_pstm_nmo_X1561'      
  else 
    echo "Some strange error" 
  fi  

filein=$segypath/$file.sgy
fileout=$supath/$file.su

segyread > $fileout tape=$filein

suwind < $fileout tmin=0 tmax=0.002 > $supath/header$afile.su

afile=`expr $afile + 1`

echo $afile

done

cd $supath
cat 3D_gathers_pstm_nmo_X*.su > 3D_gathers_pstm_nmo_Xall.su

cat header*.su > 3D_gathers_pstm_nmo_Xall_headers.su







