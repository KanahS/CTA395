#!/bin/bash -l
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -r n
#PBS -j oe
#PBS -q workq

module load python/3.10.2

cd /mnt/raid-cita/ksmith/runs/

python3 asteroseis.py
