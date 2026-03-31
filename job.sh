#!/bin/bash
#SBATCH --job-name=phi1
#SBATCH --partition=astro2_short
#SBATCH --ntasks=1                       # Single task
#SBATCH --cpus-per-task=50               # Use the maximum 40 CPUs available on one node
#SBATCH --mem=64G

THREADS=50

srun /groups/astro/tx722/Cell_division_2D/mass /groups/astro/tx722/Cell_division_2D/lyotropic_division_stress.dat -o /lustre/astro/tx722/masspy-output/sperm -t$THREADS
