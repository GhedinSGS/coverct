# coverct
<img src="https://github.com/user-attachments/assets/992605e7-10f0-4656-b902-bf7a1f31d3c6" alt="coverct_logo" width="300"> <br>
CoverCt predicts (relative) viral load in a sequencing batch with state of the art performance in a full range Ct values (cycle threshold values from diagnostic quantitative PCR). Through DuckDB SQL-like on-disk processing, it is speedy and memory efficient running within seconds on a batch of samples. It summarizes per-position per-sample f
eatures into sample level summary features, which are then fed into a random forest regression (a form of classical machine learning) to interpolate missing Ct values. The model performance is explainable and input features are based on uneveneness of coverage depths across nucleotide positions and number of minor variants (MAF < 0.50) called per sample.

The inputs are designed to be either fastq, bam or vcf. On the vcf input, it calculates the unevenness in depths at mutation positions instead of every position genome-wide. This even more lightweight analysis can be done from the NCBI Activ TRACE initiative's Google Cloud pre-calculated VCF files on any SRA-deposited SARS-CoV-2 raw sequences.
