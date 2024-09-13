#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=2:00:00          	# Limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1            		# Number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=1             		# Number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=51      	# Number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=20G               	# Specify the real memory required per node. (20G)
#SBATCH --job-name getData		# You can give your job a name for easier identification (same as -J)
#SBATCH -A egl
########## Load Modules #########

module purge
#module load GNU/6.4.0-2.28  OpenMPI/2.1.2-CUDA 
module load  SRA-Toolkit/3.0.10-gompi-2023a
module load parallel/20230722-GCCcore-12.3.0

###
donwload_paired() {
  i=$(echo $1 | cut -f1)
  echo "$i"
  fastq-dump --gzip --split-files $1;
  
}

export -f donwload_paired

# Read list of SRR and create tem list of single and paired
grep -w 'Paired' SampleList.txt | cut -f1 - > tem_PE;

# Call downloading function
parallel -j 10 donwload_paired :::: tem_PE;

rm tem_PE;
