#!/bin/sh
### Set the job name
#PBS -N ige_100_500
###Set the project name, your department dc by default
#PBS -P cse
### Request email when job begins and ends
#PBS -m bea
### Specify email address to use for notification.
#PBS -M cs1120240@iitd.ac.in
####
#PBS -l select=16:ncpus=8
### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=36:00:00
#PBS -o output_ige_100_500.txt
#PBS -e errorLog_ige_100_500.txt
#### Get environment variables from submitting shell
#PBS -V
#PBS -l software=matlab
# After job starts, must goto working directory. 
# $PBS_O_WORKDIR is the directory from where the job is fired. 
cd ~/btp/global
#job
matlab -r master_ige_simple_100_500
#NOTE
# The job line is an example : users need to change it to suit their applications
# The PBS select statement picks n nodes each having m free processors
# OpenMPI needs more options such as $PBS_NODEFILE
