#!/bin/bash

ZZX=1.0
ZZY=1.73206
REDUCFACT=1
GEO=1
MAXPTS=56000


BASEFOLDER=/home/powersr/projects/
PROJFOLDER=bubbles/newlow/
PROJNAME=E_+0.0110_B_+0.000_actabs5.conf00

# BASEFOLDER=/home/powersr/projects/gen-trans/res/
# PROJFOLDER=BUBBLES_RIBBON_1e-06/AC240_rect_lat_L_72_120_membrane_bub_R_36.6_H_6.0_1x1_xyf_0.0_rf_0.0_zf_0.0_if_0.00/abs_pot/
# PROJNAME=E_+0.0600_B_+0.000_abs5.conf00


BASENAME=$BASEFOLDER$PROJFOLDER$PROJNAME

DOSNAME=$BASENAME.ldos
#CURRENTS="l0 l2 l3 l1 l4 l5 multi"
# CURRENTS="l0 l2 l3 l1 l4 "


CURRENTS="l0"
# CURRENTS="multi"

# CURRENTS="l0"
# CURRENTS="multi"
# CURRENTS="l0 l2 l3 l1 multi"

NUMPTS=`cat $DOSNAME | awk 'END{print NR}'`
MAXDOS=`cat $DOSNAME | awk 'BEGIN{max=0} {if ($3>max) max=$3} END{printf "%lf",  max}'`
MAXX=`cat $DOSNAME | awk 'BEGIN{max=0} {if ($1>max) max=$1} END{print max}'`
MINX=`cat $DOSNAME | awk 'BEGIN{min=1000} {if ($1<min) min=$1} END{print min}'`
MAXY=`cat $DOSNAME | awk 'BEGIN{max=0} {if ($2>max) max=$2} END{print max}'`
MINY=`cat $DOSNAME | awk 'BEGIN{min=1000} {if ($2<min) min=$2} END{print min}'`
DMAX=`cat $DOSNAME | awk 'BEGIN{sum=0} {sum += $3} END{print 2*sum/NR}'`

XE=$ZZX
YE=$ZZY

if [ $GEO -eq 1 ]
then
	XE=$ZZY
	YE=$ZZX
fi

echo "$NUMPTS points..."
RESCALE=1

if [ $NUMPTS -le $MAXPTS ]
then 
    echo "OK! Very little rescaling required!"
    XCELLS=`echo "$MAXX/$XE - $MINX/$XE + 1"  | bc `
    YCELLS=`echo "$MAXY/$YE - $MINY/$YE + 1"  | bc `
    NUMPTS=$(($XCELLS*$YCELLS))
    echo "Rescaling... to "$NUMPTS" points"
fi

if [ $NUMPTS -ge $MAXPTS ]
then

  while [ $NUMPTS -ge $MAXPTS ]
  do 
	  
	  XE=`echo "$XE * $RESCALE"  | bc -l `
	  YE=`echo "$YE * $RESCALE"  | bc -l `
	  XCELLS=`echo "$MAXX/$XE - $MINX/$XE + 1"  | bc `
	  YCELLS=`echo "$MAXY/$YE - $MINY/$YE + 1"  | bc `
	  NUMPTS=$(($XCELLS*$YCELLS))
	  echo "Rescaling... to "$NUMPTS" points"
	  RESCALE=$(($RESCALE*2))
  done
fi


echo $DMAX
echo $XE $YE
echo $XCELLS $YCELLS




grep -v ^# ${DOSNAME}  | awk -vd1=$XCELLS -vd2=$YCELLS -vxe=$XE -vye=$YE -vxmin=$MINX -vymin=$MINY 'BEGIN {for(i=0; i<d1*d2; i++) arr[i]=0.0} {arr[int(  ($1 - xmin ) /xe)*d2 + int(($2 - ymin ) /ye )] += $3 } END {for(i=0; i<d1; i++) for(j=0; j<d2; j++) printf "%e %e %e\n", xmin + i*xe , ymin+j*ye,  arr[d2*i+j]} ' >  ${DOSNAME}_red.dat

# echo ${DOSNAME}_red.dat
NEWMAXDOS=`cat ${DOSNAME}_red.dat | awk 'BEGIN{max=0} {if ($3>max) max=$3} END{printf "%lf",  max}'`



# Files containing reduced current and DOS maps
for CSEL in $CURRENTS; do 
CNAME=$BASENAME.cmaps_${CSEL}
grep -v ^# ${CNAME}  | awk -vd1=$XCELLS -vd2=$YCELLS -vxe=$XE -vye=$YE -vxmin=$MINX -vymin=$MINY 'BEGIN {for(i=0; i<d1*d2; i++) arr[i]=0.0; arr2[i]=0.0} {arr[int(  ($1 - xmin ) /xe)*d2 + int(($2 - ymin ) /ye )] += $3; arr2[int(  ($1 - xmin ) /xe)*d2 + int(($2 - ymin ) /ye )] += $4 } END {for(i=0; i<d1; i++) for(j=0; j<d2; j++) printf "%e %e %e %e\n", xmin + i*xe , ymin+j*ye,  arr[d2*i+j], arr2[d2*i+j]} ' > ${CNAME}_red.dat;
done



#HALL BAR CURRENT SUMS
FILELIST=`for CSEL in $CURRENTS; do CNAME=$BASENAME.cmaps_${CSEL}; echo -n ${CNAME}_red.dat" "; done`
# echo $FILELIST
paste $FILELIST | awk '{print $1, $2, $3+$7+$11+$15+$19+$23, $4+$8+$12+$16+$20+$24 }' > ${BASENAME}_allcurrentsum.dat
paste $FILELIST | awk '{print $1, $2, $3+$7+$11, $4+$8+$12 }' > ${BASENAME}_halfcurrentsum.dat






