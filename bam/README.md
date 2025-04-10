Running preprocess_depth_var.sh executes the flu_pipeline that runs bwa mapping raw reads to a reference (such as influenza) which is based on the Anderson lab's repo: https://github.com/andersen-lab/avian-influenza. There are commented lines that are for downsampling purposes. In effect, this pipeline takes fastq inputs and outputs the samtools mpileup depth per genome position and the ivar called variants at a minimum MAF of 0.001.

These are the inputs to the coverCt model script in the model/ directory.
