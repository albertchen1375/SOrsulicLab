#!/bin/bash
#$ -cwd                      # Run in current working directory
# error = Merged with joblog
#$ -o joblog.$JOB_ID         # Output log file
#$ -j y                      # Merge error and output logs
#$ -l h_rt=48:00:00,h_data=48G  # Request 48 hours and 16GB RAM
#$ -pe shared 2             # Request 2 CPU cores
#$ -M albertchen1375@ucla.edu  # Email notifications
#$ -m bea                    # Notify on (b)egin, (e)nd, (a)bort

# Echo job info
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date`
echo " "

source /u/local/Modules/default/init/bash

# Load modules
module purge
module load R/4.0.2
module load gcc/11.3.0

# Debug: Print loaded modules
echo "Loaded modules:"
module list

# Run R script
Rscript ssGSEA2.0.R \
  --input.ds="/u/home/a/albertch/sorsulic/TCGA_RNA_Seq_unique.gct" \
  --output.prefix="/u/home/a/albertch/sorsulic/results" \
  --gene.set.databases="/u/home/a/albertch/sorsulic/msigdb.v2024.1.Hs.symbols.gmt" \
  --sample.norm.type="rank" \
  --output.score.type="NES" \
  --weight=0.75 \
  --statistic="area.under.RES" \
  --nperm=1000 \
  --combine.mode="combine.add" \
  --correl.type="z.score" \
  --min.overlap=10 \
  --par=True \
  --spare.cores=1

# Echo job completion info
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date`
echo " "
