#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --account=ctb-dionneg1
#SBATCH --time=2:00:00
#SBATCH --array=1-10

i=1

for size in 20
do
for nbstudents in 50 100
do
    for nbschools in 6
    do 
        for nbextracapa in 100
        do 
            for seed in 0 1000 2000 3000 4000
            do 
                if [ $SLURM_ARRAY_TASK_ID -eq $i ]
                then
                    python rd_instance_generator.py -size_boots $size -nbschools $nbschools -nbstudents $nbstudents -nbextracapa $nbextracapa -seed $seed
                fi
                (( i = $i +1 ))
            done 
        done
    done
done
done
