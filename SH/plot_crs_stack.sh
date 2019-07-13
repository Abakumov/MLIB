#!/bin/bash



# PLOT CRS STACKs and Attributes

# Plot Parameter
CLIP_stack=0.1
CLIP_Rnip=45000
CLIP_Rn=2000
CLIP_angle=60
CLIP_coher=0.05



labelsize=9
#titlesize=12

hbox_stack=5.1
#wbox=$((1.63 * $hbox))
wbox_stack=9.1




hbox=4.5
#wbox=$((1.63 * $hbox))
wbox=9.1

D1S=2.5		# Quality of the picture, small is high
D2S=2.5 	# Quality of the picture, small is high
axeswidth=0.5


# ZOOM Parameter
cdpmin=176
cdpmax=8126

tmin=0
tmax=6

#d2num=$((($cdpmax - $cdpmin)/20))
d2num=300
d1num=1
#d1num=$((($tmax - $tmin)/10))


# Legend parameter
legend=1
lx=$(echo "$wbox*0.5" | bc -l)
ly=$(echo "$hbox*0.25" | bc -l)
lwidth=$(echo "$wbox*0.35" | bc -l)
lheight=$(echo "$hbox*0.015" | bc -l)




# Time sections colorbar (rgb2) 
BRGB2="0.878,0.000,0.314"
GRGB2="0.941,0.941,0.941"
WRGB2="0.000,0.251,0.753"

BHLS="0.000,0.500,1.000"
GHLS="0.333,0.500,1.000"
WHLS="0.667,0.500,1.000"



# CRS FILES 
file_name=dobre.noise_muted
location_crs=/scratch/local1/u300144/CRS_TEST/OUTPUT_CRS


# PDF FILES
location_pdf=/scratch/local1/u300144/plots_crs



# Plot OPTSTACK
# -------------
input_su=$location_crs/$file_name.Optstack
output_pdf=$location_pdf/$file_name.cmps$cdpmin-$cdpmax.Optstack.pdf

segyclean < $input_su |\
supsimage \
	f2=1 d2=1 n2tic=10 \
	d2num=$d2num \
	d1s=$D1S d2s=$D2S \
	d1num=$d1num n1tic=5 \
	brgb=$BRGB2 grgb=$GRGB2 wrgb=$WRGB2 bps=24 \
	clip=$CLIP_stack wbox=$wbox_stack hbox=$hbox_stack \
	f2=$cdpmin \
	labelsize=$labelsize \
	titlesize=$titlesize \
	axeswidth=$axeswidth \
	label1="TWT [s]" label2="CMP-Nummer" \
	curve=box1_crs_zoom.txt,box2_crs_zoom.txt,box3_crs_zoom.txt npair=5,5,5 \
	curvewidth=0.7 \
> tmp.ps

#gimp tmp.ps
ps2pdf	-dPDFSETTINGS=/default -dUseFlateCompression=true \
	-dOptimize=true -dEPSCrop tmp.ps 


pdfcrop	tmp.pdf $output_pdf

#pdfedit $output_pdf

rm tmp*
acroread $output_pdf 





# Plot COHERENCY
# ---------
input_su=$location_crs/$file_name.Optcoher_1
output_pdf=$location_pdf/$file_name.cmps$cdpmin-$cdpmax.Optcoher_1.pdf
rm $output_pdf

segyclean < $input_su |\
supsimage \
	f2=1 d2=1 n2tic=10 \
	d2num=$d2num \
	d1s=$D1S d2s=$D2S \
	d1num=$d1num n1tic=5 \
	bhls=$BHLS ghls=$GHLS whls=$WHLS bps=24 \
	clip=$CLIP_coher wbox=$wbox hbox=$hbox \
	f2=$cdpmin \
	labelsize=$labelsize \
	titlesize=$titlesize \
	axeswidth=$axeswidth \
	label1="TWT [s]" label2="CMP-Nummer" \
	legend=$legend \
	lx=$lx ly=$ly units="Normalisierte Amplitude" lstyle=horibottom \
	lwidth=$lwidth lheight=$lheight \
> tmp_ps

# Create the PDF of the whole data
ps2pdf -dPDFSETTINGS=/default -dUseFlateCompression=true -dOptimize=true -dEPSCrop tmp_ps tmp.pdf
pdfcrop tmp.pdf $output_pdf
rm tmp_ps tmp.pdf
#acroread $output_pdf





# Plot Rnip
# ---------
input_su=$location_crs/$file_name.OptRnip_1
output_pdf=$location_pdf/$file_name.cmps$cdpmin-$cdpmax.OptRnip_1.pdf

segyclean < $input_su |\
supsimage \
	f2=1 d2=1 n2tic=10 \
	d2num=$d2num \
	d1s=$D1S d2s=$D2S \
	d1num=$d1num n1tic=5 \
	brgb=$BRGB2 grgb=$GRGB2 wrgb=$WRGB2 bps=24 \
	clip=$CLIP_Rnip wbox=$wbox hbox=$hbox \
	f2=$cdpmin \
	labelsize=$labelsize \
	titlesize=$titlesize \
	axeswidth=$axeswidth \
	label1="TWT [s]" label2="CMP-Nummer" \
	legend=$legend \
	lx=$lx ly=$ly units="Radius R-NIP [m]" lstyle=horibottom \
	lwidth=$lwidth lheight=$lheight \
> tmp_ps

# Create the PDF of the whole data
ps2pdf -dPDFSETTINGS=/default -dUseFlateCompression=true -dOptimize=true -dEPSCrop tmp_ps tmp.pdf
pdfcrop tmp.pdf $output_pdf
rm tmp_ps tmp.pdf
#acroread $output_pdf





# Plot Rn
# -------
input_su=$location_crs/$file_name.OptRn_1
output_pdf=$location_pdf/$file_name.cmps$cdpmin-$cdpmax.OptRn_1.pdf

segyclean < $input_su |\
supsimage \
	f2=1 d2=1 n2tic=10 \
	d2num=$d2num \
	d1s=$D1S d2s=$D2S \
	d1num=$d1num n1tic=5 \
	brgb=$BRGB2 grgb=$GRGB2 wrgb=$WRGB2 bps=24 \
	clip=$CLIP_Rn wbox=$wbox hbox=$hbox \
	f2=$cdpmin \
	labelsize=$labelsize \
	titlesize=$titlesize \
	axeswidth=$axeswidth \
	label1="TWT [s]" label2="CMP-Nummer" \
	legend=$legend \
	lx=$lx ly=$ly units="Radius R-n [m]" lstyle=horibottom \
	lwidth=$lwidth lheight=$lheight \
> tmp_ps

# Create the PDF of the whole data
ps2pdf -dPDFSETTINGS=/default -dUseFlateCompression=true -dOptimize=true -dEPSCrop tmp_ps tmp.pdf
pdfcrop tmp.pdf $output_pdf
rm tmp_ps tmp.pdf
#acroread $output_pdf





# Plot angle alpha
# ----------------
input_su=$location_crs/$file_name.Optangle_1
output_pdf=$location_pdf/$file_name.cmps$cdpmin-$cdpmax.Optangle_1.pdf

segyclean < $input_su |\
supsimage \
	f2=1 d2=1 n2tic=10 \
	d2num=$d2num \
	d1s=$D1S d2s=$D2S \
	d1num=$d1num n1tic=5 \
	brgb=$BRGB2 grgb=$GRGB2 wrgb=$WRGB2 bps=24 \
	clip=$CLIP_angle wbox=$wbox hbox=$hbox \
	f2=$cdpmin \
	labelsize=$labelsize \
	titlesize=$titlesize \
	axeswidth=$axeswidth \
	label1="TWT [s]" label2="CMP-Nummer" \
	legend=$legend \
	lx=$lx ly=$ly units="Winkel [Grad]" lstyle=horibottom \
	lwidth=$lwidth lheight=$lheight \
> tmp_ps

# Create the PDF of the whole data
ps2pdf -dPDFSETTINGS=/default -dUseFlateCompression=true -dOptimize=true -dEPSCrop tmp_ps tmp.pdf
pdfcrop tmp.pdf $output_pdf
rm tmp_ps tmp.pdf
#acroread $output_pdf

ls /scratch/local1/u300144/plots_crs/ -ltrh


exit
