#!/bin/bash
#SBATCH --partition=bgmp            ### Partition (like a queue in PBS)
#SBATCH --job-name=dedupe        ### Job Name
#SBATCH --output=output/dedupe_%j.out        ### File in which to store job output
#SBATCH --error=output/dedupe_%j.err         ### File in which to store job error messages
#SBATCH --time=0-01:00:00           ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                   ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1         ### Number of tasks to be launched per Node
#SBATCH --account=bgmp              ### Account used for job submission
#SBATCH --cpus-per-task=8

sam_in=$1
runnum=$2

#enter environment 
conda activate bgmp_py310
python --version 
#run python code
/usr/bin/time -v python li_deduper.py -f ${sam_in} -o file_deduped${runnum}.sam 
echo "deduped"
exit