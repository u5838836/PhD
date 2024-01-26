#!/bin/bash
#PBS -P vq37
#PBS -q hugemembw
#PBS -l walltime=48:00:00
#PBS -l mem=1020GB
#PBS -l wd
#PBS -l ncpus=28

module load R/3.6.1
module load intel-compiler/2020.2.254
module load intel-inspector/2019.5.0.602103
module load intel-ipp/2020.2.254
module load intel-itac/2020.2.031 
module load gcc/system
Rscript /home/138/xn2595/Spatial/bin_infill_full.R
