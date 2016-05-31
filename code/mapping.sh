#!/bin/bash

ZZX=1.0
ZZY=1.73206
REDUCFACT=1
GEO=0
MAXPTS=10000

DOSNAME=/home/spow/projects/gen_trans/res/ANTIDOTS_RIBBON_1e-06/ZZ140_trig_lat_L_7_circ_dot_R_3.0_10x4_xyf_0.0_rf_0.0/E_+0.20_B_+0.000_run.conf00.ldos
CURRENTNAME=/home/spow/projects/gen_trans/res/ANTIDOTS_RIBBON_1e-06/ZZ140_trig_lat_L_7_circ_dot_R_3.0_10x4_xyf_0.0_rf_0.0/E_+0.20_B_+0.000_run.conf00.cmaps_l0
NUMPTS=`cat $DOSNAME | awk 'END{print NR}'`
MAXDOS=`cat $DOSNAME | awk 'BEGIN{max=0} {if ($3>max) max=$3} END{print max}'`
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

# if [ $NUMPTS -le $MAXPTS ]
# then 
# 
# 
# 
# 
# 
# fi

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

# XCELLS=`echo "$MAXX/$XE - $MINX/$XE + 1"  | bc `
# YCELLS=`echo "$MAXY/$YE - $MINY/$YE + 1"  | bc `

echo $DMAX
echo $XE $YE
echo $XCELLS $YCELLS






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





 
