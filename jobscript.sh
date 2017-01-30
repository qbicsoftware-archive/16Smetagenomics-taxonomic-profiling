#!/bin/sh
#PBS -q cfc
#PBS -A qbic
#PBS -l nodes=1:ppn=10:cfc
#PBS -l walltime=40:00:00
#PBS -e ../logs/jobscript.{job.rule.name}.e$PBS_JOBID
#PBS -o ../logs/jobscript.{job.rule.name}.o$PBS_JOBID
# properties = {properties}

module load qbic/anaconda
module load devel/java_jdk/1.8.0u121   
module load qbic/clipandmerge/1.7.5
module load qbic/malt/0.3.8

{exec_job}
exit 0
