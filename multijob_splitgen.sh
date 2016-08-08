#!/bin/bash
#takes a limited number of options and assumes the rest can be fed directly to the main program as input
#this version is customised for a sublattice interface system ...is it?
#everything from 5th argument on is passed to program
#this version splits the generation and calculation processes

CFGFILE=$1 	#path of main cfgfile to be used
NUMCONFS=$2
NUMPROCS=$3	#number of processes to run per node/submission
PROCSPCONF=$4	#total number of processes/tasks per configuration
WALLT=72:00:00
GENWALLT=6:00:00	#Time allowed for generation process - probably too long...


ACTPROCSPCONF=$(($PROCSPCONF + 1))

NUMNODES=$(($PROCSPCONF/$NUMPROCS))

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
echo "# Processors per job: $NUMPROCS"
echo "# Number of jobs: $NUMNODES"   


echo ""

echo "using config file: $CFGFILE  ($FILEREDUX)"
echo ""
echo " ### start of config file"
cat $1
echo " ### end of config file"
echo ""

echo "#additional command line variables: " $CMDBITS
echo ""

CMDFILE=`cat $1 | tr '\n' ' '`
echo "COMMAND LINE: ./test" $CMDFILE $CMDBITS



mkdir -p subfiles

#loop over the various configurations
for THISCONF in `seq 0 $(($NUMCONFS -1))`; do

#first generate the scripts for proc 0 - this generates the geometries
cat > subfiles/$NAME.conf_$THISCONF.gen.bat << EOF
#!/bin/sh
#PBS -q fotonano
#PBS -l nodes=1:ppn=1
#PBS -l walltime=${GENWALLT}
#PBS -o output/\$PBS_JOBNAME.\$PBS_JOBID
#PBS -m a
#PBS -j oe

source ~/.bashrc

module load npa-cluster-setup
module load intel/2015.1.133

ulimit -s unlimited
date

cd \$PBS_O_WORKDIR

d=.src_\$PBS_JOBID
trap "rm -rf \$d" SIGINT SIGKILL SIGTERM
mkdir \$d
cd \$d
cp -rf ../code/* .
make -f cluster.mk clean
make -f cluster.mk
exe=./test 	

export OMP_NUM_THREADS=1

\$exe $CMDFILE $CMDBITS -splitgen 1 -confnum $THISCONF -numofprocs $ACTPROCSPCONF -thisproc 0

cd ../

rm -rf \$d

EOF


#submit generation (task 0) job and store ID
# JOBNUM=123456
JOBNUM=$(qsub subfiles/$NAME.conf_$THISCONF.gen.bat)


#iterator used in node loops.
TASKIT=1

#loop over each node/job submission
for THISNODE in `seq 0 $(($NUMNODES -1))`; do

cat > subfiles/$NAME.conf_$THISCONF.node_$THISNODE.bat << EOF
#!/bin/sh
#PBS -q fotonano
#PBS -l nodes=1:ppn=${NUMPROCS}
#PBS -l walltime=${WALLT}
#PBS -o output/\$PBS_JOBNAME.\$PBS_JOBID
#PBS -m a
#PBS -j oe

source ~/.bashrc

module load npa-cluster-setup
module load intel/2015.1.133

ulimit -s unlimited
date

cd \$PBS_O_WORKDIR

d=.src_\$PBS_JOBID
trap "rm -rf \$d" SIGINT SIGKILL SIGTERM
mkdir \$d
cd \$d
cp -rf ../code/* .
make -f cluster.mk clean
make -f cluster.mk
exe=./test 	

export OMP_NUM_THREADS=1

EOF


#tasks per job loop
for THISTASK in `seq 0 $(($NUMPROCS -1))`; do

cat >> subfiles/$NAME.conf_$THISCONF.node_$THISNODE.bat << EOF
\$exe $CMDFILE $CMDBITS -splitgen 1 -confnum $THISCONF -numofprocs $ACTPROCSPCONF -thisproc $TASKIT &
EOF


TASKIT=$(($TASKIT + 1))
#end of tasks loop
done

cat >> subfiles/$NAME.conf_$THISCONF.node_$THISNODE.bat << EOF
wait
cd ../
rm -rf \$d
EOF


qsub -W depend=afterok:$JOBNUM subfiles/$NAME.conf_$THISCONF.node_$THISNODE.bat
#end of node/job loop
done



#end configuration loop
done

# 
# 
# 
# 
# 
# for FILENUM in `seq 0 $(($NUMNODES -1))`; do
# cat > subfiles/$NAME.part_$FILENUM.bat << EOF
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
# cp -rf ../code/* .
# make -f cluster.mk clean
# make -f cluster.mk
# exe=./test 	
# 
# #NPROCS=\$(wc -l \$PBS_NODEFILE)
# #export OMPI_MCA_mpi_paffinity_alone=1
# export OMP_NUM_THREADS=1
# 
#
# 
# EOF
# 
# done
# 
# 
# #loop over different (disorder) configurations
# for THISCONF in `seq 0 $(($NUMCONFS -1))`; do
# 
#   #loop over different tasks (usually energy, b-field, or position ranges) within one configuration
#   for THISCONFTASK in `seq 0 $(($PROCSPCONF - 1))`; do
#   
#     THISNUMTASKS=$(($NUMPROCS * $OVERLOAD))
#     if [ $THISSUBMISSION -eq $(($NUMNODES - 1)) ]
#     then
# 	if [ $REMAINDER -gt 0 ]
# 	then
# 	THISNUMTASKS=$REMAINDER
# 	fi
#     fi
# 
# if [ $THISSUBTASK -ne $(($THISNUMTASKS - 1)) ]
# then
# cat >> subfiles/$NAME.part_$THISSUBMISSION.bat << EOF
# \$exe $CMDFILE $CMDBITS -confnum $THISCONF -numofprocs $PROCSPCONF -thisproc $THISCONFTASK &
# EOF
# fi
# 
# 
# if [ $THISSUBTASK -eq $(($THISNUMTASKS - 1)) ]
# then
# cat >> subfiles/$NAME.part_$THISSUBMISSION.bat << EOF
# \$exe $CMDFILE $CMDBITS -confnum $THISCONF -numofprocs $PROCSPCONF -thisproc $THISCONFTASK
# EOF
# fi
# 
# 
# if [ $(($THISSUBTASK % $NUMPROCS)) -eq $(($NUMPROCS-1)) ]
# then
# cat >> subfiles/$NAME.part_$THISSUBMISSION.bat << EOF
# wait
# EOF
# fi
# 
# 
#     THISTASK=$(($THISTASK + 1))
#     THISSUBTASK=$(($THISSUBTASK + 1))
# 
#     if [ $(($THISTASK % $TASKSPERNODE)) -eq 0 ]
#     then
#       THISSUBMISSION=$(($THISSUBMISSION + 1))
#       THISSUBTASK=0
#     fi
#    
#   done
# done
# 
# for FILENUM in `seq 0 $(($NUMNODES -1))`; do
# cat >> subfiles/$NAME.part_$FILENUM.bat << EOF
# wait
# cd ../
# rm -rf \$d
# EOF
# done
# 
# 
# for FILENUM in `seq 0 $(($NUMNODES -1))`; do
# qsub subfiles/$NAME.part_$FILENUM.bat
# done
# 
#    
#     
#     
#     
#     
# 
# 
# 
# 
# 
# 
# 
