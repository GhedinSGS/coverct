#!/bin/bash
# Define paths
PROJECT=$PWD
COW=$PROJECT/fastq/
OUTPUT_DIR=$PROJECT/out2/
REF=~/code/bin/flu/reference.fasta
BIN=~/code/bin

# Create swarm file with headers
echo "#SWARM --logdir swarmlogs --module samtools,bwa,java,ivar --gb-per-process 8 -t 4" > cow_ds.swarm

# Process paired-end files with _1/_2 pattern
find "$COW" -name "*_1.fastq.gz" | while read -r read1; do
    # Skip if this is part of _READ1 pattern (already processed)
    if [[ "$read1" != *"_READ1"* ]]; then
        read2="${read1/_1/_2}"
        if [ -f "$read2" ]; then
            echo "flu_pipeline.sh ${read1} ${read2} ${OUTPUT_DIR} ${REF} ${REF} ${BIN}" >> cow_ds.swarm
        else
            echo "Warning: No matching _2 file found for ${read1}"
        fi
    fi
done
