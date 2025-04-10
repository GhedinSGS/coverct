#!/bin/bash
# Check if required arguments are provided - make argument check flexible for SE/PE
if [ $# -lt 5 ] || [ $# -gt 6 ]; then
    echo "Usage for paired-end reads: $0 <read1.fastq.gz> <read2.fastq.gz> <output_dir> <reference_index> <reference_fasta> <bin_dir>"
    echo "Usage for single-end reads: $0 <reads.fastq.gz> <output_dir> <reference_index> <reference_fasta> <bin_dir>"
    echo "Example PE: $0 sample_R1.fastq.gz sample_R2.fastq.gz /path/to/output/dir /path/to/bwa/index /path/to/reference.fasta /path/to/bin"
    echo "Example SE: $0 sample.fastq.gz /path/to/output/dir /path/to/bwa/index /path/to/reference.fasta /path/to/bin"
    exit 1
fi

# Load required modules
module load samtools ivar bwa java

# Determine if input is single or paired end and set variables accordingly
if [[ $# -eq 6 ]]; then
    # Paired-end case
    READ1=$1
    READ2=$2
    OUTPUT_DIR=$3
    INDEX=$4
    REF=$5
    BIN=$6
    IS_PAIRED=true
else
    # Single-end case
    READ1=$1
    OUTPUT_DIR=$2
    INDEX=$3
    REF=$4
    BIN=$5
    IS_PAIRED=false
fi

# Create output directories
mkdir -p "$OUTPUT_DIR"/{cns_seqs,depth,variants,bam,ds}
mkdir -p "$OUTPUT_DIR/ds"/{cns_seqs,depth,variants}

# Get sample name from filename - adapted for both SE and PE
n=$(basename "$READ1")
if [[ "$IS_PAIRED" == true ]]; then
    n=${n/_READ1*/}
    n=${n/_R1*/}
    n=${n/_1.fastq*/}
    n=${n/_L001*/}
else
    n=${n/.fastq*/}
    n=${n/_L001*/}
fi

echo "Processing sample: $n"
echo "Read1: $READ1"
if [[ "$IS_PAIRED" == true ]]; then
    echo "Read2: $READ2"
    echo "Running in paired-end mode"
else
    echo "Running in single-end mode"
fi


# Define BAM files
BAM="$OUTPUT_DIR/bam/${n}.bam"
DEDUP_BAM="$OUTPUT_DIR/bam/${n}.dedup.bam"
# READS1=20000
# READS2=10000
# DS1_BAM="$OUTPUT_DIR/ds/${n}-${READS1}.bam"
# DS2_BAM="$OUTPUT_DIR/ds/${n}-${READS2}.bam"

# Generate original BAM if it doesn't exist
if [ ! -f "$BAM" ]; then
    echo "Generating original BAM..."
    if [[ "$IS_PAIRED" == true ]]; then
        bwa mem -t 4 "$INDEX" "$READ1" "$READ2" | \
        #samtools view -F 4 -b | \
        samtools collate -@ 4 -O - | \
        samtools fixmate -@ 4 -m - - | \
        samtools sort -@ 4 -m 4G - -o $BAM -O BAM
    else
        bwa mem -t 4 "$INDEX" "$READ1" | \
        #samtools view -F 4 -b | \
        samtools collate -@ 4 -O - | \
        samtools fixmate -@ 4 -m - - | \
        samtools sort -@ 4 -m 4G - -o $BAM -O BAM
    fi
    samtools index "$BAM"
else
    echo "Original BAM exists, skipping alignment..."
fi

# Generate deduplicated BAM if it doesn't exist
if [ ! -f "$DEDUP_BAM" ]; then
    echo "Generating deduplicated BAM..."
    samtools markdup -@ 4 -r -l 20000 $BAM "$DEDUP_BAM"
    samtools index "$DEDUP_BAM"
else
    echo "Deduplicated BAM exists, skipping deduplication..."
fi

# # Generate downsampled BAMs if they don't exist
# if [ ! -f "$DS1_BAM" ] || [ ! -f "$DS2_BAM" ]; then
#     echo "Performing downsampling..."
#     # Calculate fractions for downsampling
#     fraction1=$(samtools idxstats "$BAM" | cut -f3 | awk -v ct=$READS1 'BEGIN {total=0} {total += $1} END {print ct/total}')
#     fraction2=$(awk -v fraction1="$fraction1" -v reads1=$READS1 -v reads2=$READS2 'BEGIN {print fraction1 * (reads2/reads1)}')

#     # Create first downsample if it doesn't exist
#     if [ ! -f "$DS1_BAM" ]; then
#         java -jar "$BIN/picard.jar" DownsampleSam P=$fraction1 I="$BAM" O="$DS1_BAM"
#         samtools index "$DS1_BAM"
#     fi

#     # Create second downsample if it doesn't exist
#     if [ ! -f "$DS2_BAM" ]; then
#         java -jar "$BIN/picard.jar" DownsampleSam P=$fraction2 I="$BAM" O="$DS2_BAM"
#         samtools index "$DS2_BAM"
#     fi
# else
#     echo "Downsampled BAMs exist, skipping downsampling..."
# fi

# Function to process a BAM file
process_bam() {
    local bam_file=$1
    local output_prefix=$2
    local output_dir=$3
    
    # Generate consensus if it doesn't exist
    if [ ! -f "${output_dir}/cns_seqs/${output_prefix}_cns.fa" ]; then
        echo "Generating consensus for ${output_prefix}..."
        samtools mpileup -A -d 0 -aa -Q 0 "$bam_file" | \
            ivar consensus -p "${output_dir}/cns_seqs/${output_prefix}_cns" -t 0.5 -m 1
    fi

    # Generate depth file if it doesn't exist
    if [ ! -f "${output_dir}/depth/${output_prefix}_depth.tsv" ]; then
        echo "Calculating depth for ${output_prefix}..."
        samtools depth -aa "$bam_file" > "${output_dir}/depth/${output_prefix}_depth.tsv"
    fi

    # Call variants if they don't exist
    if [ ! -f "${output_dir}/variants/${output_prefix}_variants.tsv" ]; then
        echo "Calling variants for ${output_prefix}..."
        samtools mpileup --reference "$REF" -A -d 0 -aa -Q 0 "$bam_file" | \
            ivar variants -p "${output_dir}/variants/${output_prefix}_variants" -t 0.001 -m 1
    fi
}

# Process original BAM
echo "Processing original BAM..."
process_bam "$BAM" "${n}" "$OUTPUT_DIR"

# Process deduplicated BAM
echo "Processing deduplicated BAM..."
process_bam "$DEDUP_BAM" "${n}-dedup" "$OUTPUT_DIR"

# # Process first downsample
# echo "Processing ${READS1} downsample..."
# process_bam "$DS1_BAM" "${n}-${READS1}" "$OUTPUT_DIR"

# # Process second downsample
# echo "Processing ${READS2} downsample..."
# process_bam "$DS2_BAM" "${n}-${READS2}" "$OUTPUT_DIR"

echo "Completed processing sample: $n"
echo "-----------------------------------"
echo "Sample processed successfully!"
