#!/bin/bash

# ask for 1Gb of RAM with an upper bound of 2G if you exceed h_vmem the task will be cancelled without warning
#$ -l tmem=4G
#$ -l h_vmem=4G

# ask for 1G of tmp/scratch space
#$ -l tscratch=1G

# force it to only schedule a task on machines that actually have 1Gb of tmp  free at runtime
#$ -l scratch0free=1G

# set max runtime for each task
#$ -l h_rt=48:0:0

# Merge stdout and stderr to a single output file
# you can send them to separate files and 
#$ -j y

# give the array job a name
#$ -N pfam_nw

#set number of array tasks, each task gets a value from this range 
# in environment variable SGE_TASK_ID
#$ -t 100-504

#pass your whole environment.
#$ -V

# setting a GPU and selection specfic hosts
#$ -l hostname=(zeus1.local, zeus2.local)
#$ -l gpu=true
#$ -pe gpu 1

#### Run the application.

# we setup a cleanup function in case we or the scheduler needs to kill
# our task and we need to tidy up the /tmp space we were using before
# the script terminates
cleanup() {
    # annouce what happened to stdout
    echo "Received kill signal. Cleaning up..."

    # now delete the tmp dir contents that the task was using
    cd ~/
    rm -rf /scratch0/pfam_nw_${SGE_TASK_ID}
    
    exit 0
}
# now register which signals will cause cleanup() to be invoked
trap cleanup SIGINT
trap cleanup SIGTERM
trap cleanup EXIT
trap cleanup ERR

# print the location that task is running to stdout (handy for debugging)
hostname 

# Create a directory in the tmp space that is uniquely named for this
# task using the task ID. Stops any multi-thread collisions
# though I've not done it here you could make a unique name with
# ${JOB_NAME}_${SGE_TASK_ID}. Or if you're paranoid someone will
# be using the same job name as you ${JOB_ID}_${SGE_TASK_ID}
# or you want your tmp space to be thread safe between array submissions
mkdir /scratch0/pfam_nw_${SGE_TASK_ID}/
cd /scratch0/pfam_nw_${SGE_TASK_ID}/

# Copy any files the task needs locally to the tmp space
cp /home/dbuchan/pfam_nw/random_pfam_reps.fa /scratch0/pfam_nw_${SGE_TASK_ID}/
cp /home/dbuchan/pfam_nw/${SGE_TASK_ID}_pfam_random /scratch0/pfam_nw_${SGE_TASK_ID}/

##### NOW DO STUFF
# Send the command to STDOUT as a string (handy for debugging)
# I'm only really doing one number crunching thing here but you might
# have lots to do
echo "python /home/dbuchan/profile_drift_new/scripts/rep_distance_matrix/pfam_reps_nw.py /scratch0/pfam_nw_${SGE_TASK_ID}/random_pfam_reps.fa /scratch0/pfam_nw_${SGE_TASK_ID}/${S
GE_TASK_ID}_pfam_random > /scratch0/cp /home/dbuchan/pfam_nw/random_pfam_reps.fa /scratch0/pfam_nw_${SGE_TASK_ID}/pfam_nw_${SGE_TASK_ID}/${SGE_TASK_ID}_hits.csv 2> /scratch0/pfam_nw_${SGE_TASK_ID}/${SGE_TASK_ID}_hits.err"
python /home/dbuchan/profile_drift_new/scripts/rep_distance_matrix/pfam_reps_nw.py /scratch0/pfam_nw_${SGE_TASK_ID}/random_pfam_reps.fa /scratch0/pfam_nw_${SGE_TASK_ID}/${SGE_TAS
K_ID}_pfam_random > /scratch0/pfam_nw_${SGE_TASK_ID}/${SGE_TASK_ID}_hits.csv 2> /scratch0/pfam_nw_${SGE_TASK_ID}/${SGE_TASK_ID}_hits.err

# Now I copy the results files I want to keep back to my home directory on
# morecambe
cp /scratch0/pfam_nw_${SGE_TASK_ID}/*.csv /home/dbuchan/pfam_nw/
cp /scratch0/pfam_nw_${SGE_TASK_ID}/*.err /home/dbuchan/pfam_nw/

# Now clean up your tmp dir so that you play nice with other users
# and don't fill up the node's filesystems. We could infact just call
# a variant of the cleanup() function that was defined earlier withouy
# a "I've been killed" alert
cd ~/
rm -rf /scratch0/pfam_nw_${SGE_TASK_ID}