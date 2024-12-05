#!/bin/bash
ls /data/SGSlab/tung/240901-tonkin_rerun1/pipeline_output/minimap2/*.ivar.bam >| list_bam_to_ds.txt
ls /data/SGSlab/tung/240901-tonkin_rerun2/pipeline_output/minimap2/*.ivar.bam >> list_bam_to_ds.txt
BIN=/data/SGSlab/tung/bin/
# Remove the existing ds.swarm file to start fresh
rm -f ds.swarm

# Process each line in list_bam_to_ds.txt
while read -r file; do
    # Define variables
    BAM=$file
    name=$(basename $file .ivar.bam)
    READS1=200000

    # Append commands to ds.swarm
    echo "module load samtools java; samtools index $BAM; \\
    fraction1=\$(samtools idxstats $BAM | cut -f3 | awk -v ct=$READS1 'BEGIN {total=0} {total += \$1} END {print ct/total}'); \\
    fraction2=\$(awk -v fraction1="\$fraction1" 'BEGIN {print fraction1 / 10}'); \\
    java -jar $BIN/picard.jar DownsampleSam P=\$fraction1 I=$BAM O=ds/${name}-${READS1}.bam; samtools index ds/${name}-${READS1}.bam; \\
    java -jar $BIN/picard.jar DownsampleSam P=\$fraction2 I=$BAM O=ds/${name}-${READS1}.bam; samtools index ds/${name}-${READS1}.bam;" >> ds.swarm
done < list_bam_to_ds.txt
