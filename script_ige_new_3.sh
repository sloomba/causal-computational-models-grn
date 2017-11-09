#!/bin/sh
### Set the job name
#PBS -N ige_new_3
###Set the project name, your department dc by default
#PBS -P cse
### Request email when job begins and ends
#PBS -m bea
### Specify email address to use for notification.
#PBS -M cs1120114@iitd.ac.in
####
#PBS -l select=4:ncpus=4
### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=48:00:00
#PBS -o output_ige_new_3.txt
#PBS -e errorLog_ige_new_3.txt
#### Get environment variables from submitting shell
#PBS -V
#PBS -l software=matlab
# After job starts, must goto working directory. 
# $PBS_O_WORKDIR is the directory from where the job is fired. 
cd ~/btp/global
#job
matlab -r master_ige_simple_new_3
#NOTE
# The job line is an example : users need to change it to suit their applications
# The PBS select statement picks n nodes each having m free processors
# OpenMPI needs more options such as $PBS_NODEFILE
