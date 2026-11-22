#!/bin/bash
#SBATCH --job-name=mertk-mask-reach
#SBATCH --account=pi_pallard
#SBATCH --partition=campus
#SBATCH --qos=campus24
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
set -eu
export OPENBLAS_NUM_THREADS=1
project=/u/scratch/s/sujit009/c288-mertk-binder
"$project/envs/project/bin/python" mask_geometry.py --native target.pdb --context *_full_context.pdb --binder *_binder.pdb --out mask_geometry.tsv
