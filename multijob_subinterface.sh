#!/bin/bash
#takes a limited number of options and assumes the rest can be fed directly to the main program as input
#this version is customised for a sublattice interface system ...is it?
#everything from 5th argument on is passed to program

CFGFILE=$1 	#path of main cfgfile to be used
NUMCONFS=$2
NUMPROCS=$3
PROCSPCONF=$4
OVERLOAD=1
WALLT=72:00:00

TASKS=$(($NUMCONFS * $PROCSPCONF))
NUMNODES=$(($TASKS/($NUMPROCS*$OVERLOAD)))
REMAINDER=$(($TASKS % ($NUMPROCS*$OVERLOAD)))
TASKSPERNODE=$(($NUMPROCS*$OVERLOAD))

CMDBITS=""
DATETIME=`date +%S_%N`
FILEREDUX=`echo $CFGFILE | tr  / _ | tr  . _`
NAME=$FILEREDUX"_"$DATETIME

if [ $# -ge 5 ]
then
for i in ${@:5}; do
	CMDBITS=$CMDBITS" "$i
done
fi

if [ $REMAINDER -ne 0 ]
then
  NUMNODES=$(($NUMNODES +1))
fi

echo "# configurations: $NUMCONFS"
echo "# processors per configuration: $PROCSPCONF"
echo "# tasks: $TASKS"
echo "# Processors per node: $NUMPROCS"
echo "# tasks per processor overloading: $OVERLOAD"
echo "# tasks per node: $TASKSPERNODE"
echo "# single node jobs: $NUMNODES"   

echo "REMAINDER: $REMAINDER "   

echo ""

echo "using config file: $CFGFILE  ($FILEREDUX)"
echo ""
echo " ### start of config file"
cat $1
echo " ### end of config file"
echo ""

echo "#additional command line variables: " $CMDBITS
echo ""

echo "COMMAND LINE: ./test" `cat $1` $CMDBITS

THISTASK=0
THISSUBMISSION=0
THISSUBTASK=0

mkdir -p subfiles

for FILENUM in `seq 0 $(($NUMNODES -1))`; do
cat > subfiles/$NAME.part_$FILENUM.bat << EOF
#!/bin/sh
#PBS -q fotonano
#PBS -l nodes=1:ppn=${NUMPROCS}
#PBS -l walltime=${WALLT}
#PBS -o output/\$PBS_JOBNAME.\$PBS_JOBID
#PBS -m a
#PBS -j oe

#source /zhome/0e/2/36189/.bashrc
source ~/.bashrc

#module purge
module load npa-cluster-setup
module load intel/2015.1.133

ulimit -s unlimited
date

#env
#echo ""
#env | grep NPA_


cd \$PBS_O_WORKDIR

d=.src_\$PBS_JOBID
trap "rm -rf \$d" SIGINT SIGKILL SIGTERM
mkdir \$d
cd \$d
cp -rf ../code/* .
make clean
make -f cluster.mk
cd ../
exe=\$d/test 	

#NPROCS=\$(wc -l \$PBS_NODEFILE)
#export OMPI_MCA_mpi_paffinity_alone=1
export OMP_NUM_THREADS=1

EOF

done


#loop over different (disorder) configurations
for THISCONF in `seq 0 $(($NUMCONFS -1))`; do

  #loop over different tasks (usually energy, b-field, or position ranges) within one configuration
  for THISCONFTASK in `seq 0 $(($PROCSPCONF - 1))`; do
  
    THISNUMTASKS=$(($NUMPROCS * $OVERLOAD))
    if [ $THISSUBMISSION -eq $(($NUMNODES - 1)) ]
    then
	if [ $REMAINDER -gt 0 ]
	then
	THISNUMTASKS=$REMAINDER
	fi
    fi

if [ $THISSUBTASK -ne $(($THISNUMTASKS - 1)) ]
then
cat >> subfiles/$NAME.part_$THISSUBMISSION.bat << EOF
\$exe `cat $1` $CMDBITS "-confnum $THISCONF -numofprocs $PROCSPCONF -thisproc $THISCONFTASK" &
EOF
fi


if [ $THISSUBTASK -eq $(($THISNUMTASKS - 1)) ]
then
cat >> subfiles/$NAME.part_$THISSUBMISSION.bat << EOF
\$exe `cat $1` $CMDBITS "-confnum $THISCONF -numofprocs $PROCSPCONF -thisproc $THISCONFTASK" 
EOF
fi


if [ $(($THISSUBTASK % $NUMPROCS)) -eq $(($NUMPROCS-1)) ]
then
cat >> subfiles/$NAME.part_$THISSUBMISSION.bat << EOF
wait
EOF
fi


    THISTASK=$(($THISTASK + 1))
    THISSUBTASK=$(($THISSUBTASK + 1))

    if [ $(($THISTASK % $TASKSPERNODE)) -eq 0 ]
    then
      THISSUBMISSION=$(($THISSUBMISSION + 1))
      THISSUBTASK=0
    fi
   
  done
done

for FILENUM in `seq 0 $(($NUMNODES -1))`; do
cat >> subfiles/$NAME.part_$FILENUM.bat << EOF
wait
rm -rf \$d
EOF
done


for FILENUM in `seq 0 $(($NUMNODES -1))`; do
qsub subfiles/$NAME.part_$FILENUM.bat
done

   
    
    
    
    







