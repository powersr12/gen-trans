#!/bin/bash

ZZX=1.0
ZZY=1.73206
REDUCFACT=1
GEO=0
MAXPTS=100000

XMIN=0
XMAX=850
YMIN=0
YMAX=420

BASEFOLDER=/home/ICN2/spower/DTU_GEN_RES/
PROJFOLDER=ANTIDOTS_RIBBON_1e-06/ZZ480_rect_lat_L_60_104_circ_dot_R_20.0_4x8_xyf_0.0_rf_0.0/
PROJNAME=E_+0.40_B_+41.000_NEW.conf00



BASENAME=$BASEFOLDER$PROJFOLDER$PROJNAME


DOSNAME=$BASENAME.ldos
# CURRENTS="l0 l2 l3 l1 l4 l5 multi"
CURRENTS="l0 l1"


#HALL BAR CURRENT SUMS
FILELIST=`for CSEL in $CURRENTS; do CNAME=$BASENAME.cmaps_${CSEL}; echo -n ${CNAME}_red.dat" "; done`

for FILE in $FILELIST; do 
# echo $FILE; 
cat $FILE | awk -vxmin=$XMIN -vymin=$YMIN -vxmax=$XMAX -vymax=$YMAX '{if( $1>xmin && $1 <xmax && $2>ymin && $2<ymax) print $0}' > ${FILE}_zoom.dat

done

DOSNAME=${BASENAME}.ldos_red.dat
cat ${DOSNAME}  | awk -vxmin=$XMIN -vymin=$YMIN -vxmax=$XMAX -vymax=$YMAX '{if( $1>xmin && $1 <xmax && $2>ymin && $2<ymax) print $0}' > ${BASENAME}.ldos_red.dat_zoom.dat





