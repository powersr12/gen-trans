#!/bin/bash

ZZX=1.0
ZZY=1.73206
REDUCFACT=1
GEO=0
MAXPTS=100000

# BASEFOLDER=/home/ICN2/spower/projects/gen-trans/res/
# PROJFOLDER=ANTIDOTS_RIBBON_1e-06/ZZ288_rect_lat_L_36_64_circ_dot_R_8.0_4x8_xyf_0.0_rf_0.0/
# PROJNAME=E_+0.40_B_+163.500_longbuf.conf00

BASEFOLDER=/home/ICN2/spower/projects/gen-trans/res/
PROJFOLDER=SUBLMOIRE_RIBBON_1e-06/AC110_L2_108_default_off_0.0_0.0/
PROJNAME=E_+0.0005_B_+0.000_hidef.conf00

# BASEFOLDER=/home/ICN2/spower/DTU_GEN_RES/
# PROJFOLDER=ANTIDOTS_RIBBON_1e-05/ZZ720_rect_lat_L_60_104_circ_dot_R_20.0_6x8_xyf_0.0_rf_0.0/
# PROJNAME=E_+0.40_B_+44.000_NEWlong.conf00

# BASEFOLDER=/home/spow//GEN_TRANS_RESULTS/
# PROJFOLDER=ANTIDOTS_HALLBAR_2_2_rw_6_1e-06/ZZ480_rect_lat_L_60_104_circ_dot_R_20.0_4x8_xyf_0.0_rf_0.0/
# PROJNAME=E_+0.40_B_+90.000_NEW.conf00

BASENAME=$BASEFOLDER$PROJFOLDER$PROJNAME

DOSNAME=$BASENAME.ldos
# CURRENTS="l0 l2 l3 l1 l4 l5 multi"
# CURRENTS="l0 l2 l3 l1 l4 "

 CURRENTS="l0"
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

# 
# #Create Mathematica script to generate DOS only map
# cat > tempmath.m << EOF
# 
# DOSdata =   Import["${DOSNAME}_red.dat", "Table"];
# 
# dmax = $NEWMAXDOS;
# 
# 
# DOSplot = ListDensityPlot[DOSdata, PlotRange -> {{$MINX, $MAXX}, {$MINY, $MAXY}, {0.00, dmax}}, AspectRatio -> ( $MAXY - $MINY )/( $MAXX - $MINX ), ColorFunction -> (Lighter[Blue, 1 - (#/dmax)] &), 
#   ColorFunctionScaling -> False, ClippingStyle -> {White, Blue},   MaxPlotPoints -> 200, InterpolationOrder -> 2];
# 
# Export["${DOSNAME}_map.png", DOSplot, ImageResolution -> 200];
# 
# 
# EOF


# Files containing reduced current and DOS maps
for CSEL in $CURRENTS; do 
CNAME=$BASENAME.cmaps_${CSEL}
grep -v ^# ${CNAME}  | awk -vd1=$XCELLS -vd2=$YCELLS -vxe=$XE -vye=$YE -vxmin=$MINX -vymin=$MINY 'BEGIN {for(i=0; i<d1*d2; i++) arr[i]=0.0; arr2[i]=0.0} {arr[int(  ($1 - xmin ) /xe)*d2 + int(($2 - ymin ) /ye )] += $3; arr2[int(  ($1 - xmin ) /xe)*d2 + int(($2 - ymin ) /ye )] += $4 } END {for(i=0; i<d1; i++) for(j=0; j<d2; j++) printf "%e %e %e %e\n", xmin + i*xe , ymin+j*ye,  arr[d2*i+j], arr2[d2*i+j]} ' > ${CNAME}_red.dat;
done

# #Add to Mathematica script to generate current maps with and without DOS background
# 
# cat >> tempmath.m << EOF
# 
# CurData =   Import["${CNAME}_red.dat", "Table"];
# 
# CurPlot${CSEL} =  ListVectorPlot[CurData, PlotRange -> {{$MINX, $MAXX}, {$MINY, $MAXY}}, VectorStyle -> RGBColor[0.9, 0.45, 0], VectorPoints -> {30, 30},  VectorScale -> {Medium, Scaled[1.15], Automatic}, AspectRatio -> ( $MAXY - $MINY )/( $MAXX - $MINX )];
# 
# Export["${CNAME}_map.png", CurPlot, ImageResolution -> 200];
# 
# CombinedPlot = Show[DOSplot, CurPlot];
# 
# Export["${CNAME}_map_comb.png", CombinedPlot, ImageResolution -> 200];
# 
# EOF
# 
# 
# 
# done
# 
# 
# 
# 
# cat >> tempmath.m << EOF
# 
# CombinedPlot = Show[DOSplot, CurPlotl0, CurPlotl2, CurPlotl3];
# 
# Export["${CNAME}_map_comb.png", CombinedPlot, ImageResolution -> 200];
# 
# EOF
# 


#HALL BAR CURRENT SUMS
FILELIST=`for CSEL in $CURRENTS; do CNAME=$BASENAME.cmaps_${CSEL}; echo -n ${CNAME}_red.dat" "; done`
# echo $FILELIST
paste $FILELIST | awk '{print $1, $2, $3+$7+$11+$15+$19+$23, $4+$8+$12+$16+$20+$24 }' > ${BASENAME}_allcurrentsum.dat
paste $FILELIST | awk '{print $1, $2, $3+$7+$11, $4+$8+$12 }' > ${BASENAME}_halfcurrentsum.dat











# cat > tempmath.m << EOF
# #!/bin/sh
# #PBS -q fotonano
# #PBS -l nodes=1:ppn=${NUMPROCS}
# #PBS -l walltime=${WALLT}
# #PBS -o output/\$PBS_JOBNAME.\$PBS_JOBID
# #PBS -m a
# #PBS -j oe
# 
# #source /zhome/0e/2/36189/.bashrc
# source ~/.bashrc
# 
# #module purge
# module load npa-cluster-setup
# module load intel/2015.1.133
# 
# ulimit -s unlimited
# date
# 
# #env
# #echo ""
# #env | grep NPA_
# 
# 
# cd \$PBS_O_WORKDIR
# 
# d=.src_\$PBS_JOBID
# trap "rm -rf \$d" SIGINT SIGKILL SIGTERM
# mkdir \$d
# cd \$d
# cp -rf ../src/* .
# make clean
# make
# cd ../
# exe=\$d/template 	
# 
# #NPROCS=\$(wc -l \$PBS_NODEFILE)
# #export OMPI_MCA_mpi_paffinity_alone=1
# export OMP_NUM_THREADS=1
# 
# EOF
# 





 
