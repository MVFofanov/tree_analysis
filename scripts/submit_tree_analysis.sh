#!/bin/bash

#SBATCH --job-name=snakemake_tree_analysis
#SBATCH --output=/home/zo49sog/crassvirales/phylomes/tree_analysis/slurm_logs/result_%x.%j.txt
#SBATCH --time=3-00:00:00  # This will be overwritten by the script
#SBATCH --partition=standard  # This will be overwritten by the script
#SBATCH --nodes=1  # This will be overwritten by the script
#SBATCH --ntasks=1  # This will be overwritten by the script
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB

date; hostname; pwd

# Load Snakemake environment or module
source /home/zo49sog/mambaforge/etc/profile.d/conda.sh && conda activate tree_analysis

working_dir="/home/zo49sog/crassvirales/phylomes/tree_analysis"

snakefile="/home/zo49sog/crassvirales/phylomes/tree_analysis/Snakefile"
jobs=5000
latency_wait="60"
slurm_output_dir="${working_dir}/slurm_logs/logs"
slurm_time="30:00"
slurm_partition="short"
slurm_nodes=1
slurm_ntasks=1
slurm_mem_gb=5

# Change to the working directory
cd ${working_dir}

mkdir -p ${slurm_output_dir}

snakemake --snakefile "${snakefile}" --jobs ${jobs} --cluster "sbatch --output=${slurm_output_dir}/%x_%j.out.txt --time=${slurm_time} --partition=${slurm_partition} --nodes=${slurm_nodes} --ntasks=${slurm_ntasks} --cpus-per-task={threads} --mem=${slurm_mem_gb}GB" --latency-wait ${latency_wait}

# Run Snakemake with SLURM cluster submission
#snakemake --snakefile "${snakefile}" --jobs ${jobs} --cluster "sbatch --output=${slurm_output_dir}/{rule}_%x_%j.out --time=${slurm_time} --partition=${slurm_partition} --nodes=${slurm_nodes} --ntasks=${slurm_ntasks} --cpus-per-task={threads} --mem={resources.mem_mb}MB" --latency-wait ${latency_wait}

date