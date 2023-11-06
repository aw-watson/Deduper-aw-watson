#!/bin/bash
#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp               #REQUIRED: which partition to use
#SBATCH --cpus-per-task=1                 #optional: number of cpus, default is 1
#SBATCH --mem=8GB                        #optional: amount of memory, default is 4GB
#SBATCH --output=/projects/bgmp/apwat/bioinfo/Bi624/Deduper-aw-watson/logs/deduping_choice_out%j.log
#SBATCH --error=/projects/bgmp/apwat/bioinfo/Bi624/Deduper-aw-watson/logs/deduping_choice_err%j.log
conda activate deduper

START_FILE="/projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam"
DIR="/projects/bgmp/apwat/bioinfo/Bi624/Deduper-aw-watson"

/usr/bin/time -v $DIR/watson_deduper_setup.py \
    -f $START_FILE \
    -o $DIR/input/C1_SE_uniqAlign_adjusted.sam

samtools sort -o $DIR/input/C1_SE_uniqAlign_sorted.sam $DIR/input/C1_SE_uniqAlign_adjusted.sam

/usr/bin/time -v $DIR/watson_deduper.py \
    -u STL96.txt \
    -f $DIR/input/C1_SE_uniqAlign_sorted.sam \
    -o $DIR/output/C1_SE_uniqAlign_deduped.sam \
    -m first