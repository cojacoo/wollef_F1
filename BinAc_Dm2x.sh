#!/bin/bash
#PBS -l nodes=1:ppn=3
#PBS -l walltime=72:00:00
#PBS -l mem=36gb
#PBS -S /bin/bash
#PBS -N echoRD_wollef
#PBS -j oe
#PBS -o LOG
#PBS -n

#cd /beegfs/work/ka_oj4748

echo "my Username is:"
whoami
echo "My job is running on node:"
uname -a

 
module load numlib/numpy
module load lib/matplotlib
module load lib/pandas
module load numlib/scipy

python wollef_Dm2.py &
python wollef_Em2.py &
python wollef_Fm2.py &