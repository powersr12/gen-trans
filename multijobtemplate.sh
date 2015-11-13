#!/bin/bash

LENGTH1=$1
LENGTH2=$2
GEO=$3
BUFFER=3
POSDEP=$4
DISRUNS=$5
DISPROB=$6
EDGECUT=5
SUBSTHICK=${10}
MAPS=0
EMIN=0.0
EMAX=6.0
EPTS=121
NUMCONFS=$7

NUMPROCS=${8}
PROCSPCONF=${9}
OVERLOAD=1
WALLT=72:00:00

TASKS=$(($NUMCONFS * $PROCSPCONF))
NUMNODES=$(($TASKS/($NUMPROCS*$OVERLOAD)))
REMAINDER=$(($TASKS % ($NUMPROCS*$OVERLOAD)))
TASKSPERNODE=$(($NUMPROCS*$OVERLOAD))

if [ $GEO -eq 0 ]
then
	TYPE=ZZ
fi

if [ $GEO -eq 1 ]
then
        TYPE=AC
fi
	
if [ $GEO -eq 2 ]
then
        TYPE=AC3N
fi


if [ $POSDEP -eq 0 ]
then
        PTYPE=NOPOT
fi

if [ $POSDEP -eq 1 ]
then
        PTYPE=POT
fi

if [ $POSDEP -eq 2 ]
then
        PTYPE=EPOT
	EDGECUT=15
fi



NAME=$TYPE"-"$LENGTH1"x"$LENGTH2"-"$PTYPE"-DIS"$DISRUNS"x"$DISPROB"-DI"$SUBSTHICK

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


#for i in `seq 0 $((${PROCS}-1))`; do

#echo template $1 $2 $3 $4 $3 $5 $6 $7 $8 $9 ${10} ${11} &
#done

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
cp -rf ../src/* .
make clean
make
cd ../
exe=\$d/template 	

#NPROCS=\$(wc -l \$PBS_NODEFILE)
#export OMPI_MCA_mpi_paffinity_alone=1
export OMP_NUM_THREADS=1

EOF

done



for THISCONF in `seq 0 $(($NUMCONFS -1))`; do

  

  for THISCONFTASK in `seq 0 $(($PROCSPCONF - 1))`; do
  

    THISNUMTASKS=$(($NUMPROCS * $OVERLOAD))
    if [ $THISSUBMISSION -eq $(($NUMNODES - 1)) ]
    then
        if [ $REMAINDER -gt 0 ]
        then
        THISNUMTASKS=$REMAINDER
        fi
    fi

  
   # echo $THISSUBMISSION $THISSUBTASK $THISNUMTASKS template $LENGTH $RADIUS $TBCELLS $MIDCELLS $TBCELLS $ROWS $XYERR $RERR $RDIS $PDIS $EMIN $EMAX $EPTS $THISCONF $PROCSPCONF $THISCONFTASK &

if [ $THISSUBTASK -ne $(($THISNUMTASKS - 1)) ]
then
cat >> subfiles/$NAME.part_$THISSUBMISSION.bat << EOF
\$exe $NAME $GEO $LENGTH1 $LENGTH2 $BUFFER $POSDEP $EDGECUT $SUBSTHICK $DISRUNS $DISPROB $EMIN $EMAX $EPTS $MAPS $THISCONF $PROCSPCONF $THISCONFTASK &
EOF
fi

if [ $THISSUBTASK -eq $(($THISNUMTASKS - 1)) ]
then
cat >> subfiles/$NAME.part_$THISSUBMISSION.bat << EOF
\$exe $NAME $GEO $LENGTH1 $LENGTH2 $BUFFER $POSDEP $EDGECUT $SUBSTHICK $DISRUNS $DISPROB $EMIN $EMAX $EPTS $MAPS $THISCONF $PROCSPCONF $THISCONFTASK
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
  
