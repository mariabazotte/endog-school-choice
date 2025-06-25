#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --time=15:00:00
#SBATCH --array=1-2160
#SBATCH --output=../slurm/strat-large/perf_%A_%a.out

i=1
seed_sce=0
verbose=1
nbthreads=1
timelimit=28800

nbproblems=1
nbscenarios=100 
nbvalidatescenarios=5000

for seed in 0 1000 2000 3000 4000
do
    for type in 1
    do
    for transform in 1 3
    do
    for solver in 5 7
    do
        for schools in 20 40 60
        do
        for students in 500 1000
        do 
        for seats in 1 5 15 30 45 60
        do
        for rank in 2 3 4
        do
            for warmstart in 1 
            do
                if [ $SLURM_ARRAY_TASK_ID -eq $i ]
                then
                    ./exe -type $type -transform $transform -solver $solver -instancefile "../data/instances/20bts_${schools}schools${students}students100extracapa${seed}seed.txt" -nbranks $rank -nbseats $seats -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -seed $seed_sce -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidatescenariosSAA $nbvalidatescenarios -warmstart $warmstart
                fi
                (( i = $i +1 ))
            done
        done
        done
        done
        done
    done
    done
    done
done
