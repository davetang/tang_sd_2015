# GATK pipeline

This Markdown document explains each step of the GATK pipeline that was used to process our exome data. The workflow, which follows the [GATK best practices](https://www.broadinstitute.org/gatk/guide/best-practices), has been implemented as a [Bpipe](http://docs.bpipe.org/) pipeline, where different parts of the workflow are modularised as code blocks written in the [Groovy programming language](http://www.groovy-lang.org/). A Groovy file that contains the GATK commands is used by Bpipe to run the pipeline; this document explains the GATK pipeline in the context of the Groovy file. The file begins with a title definition and the base directory is defined; ```BASEDIRNAME``` will be interpolated into the top directory when creating the Groovy file.

~~~~{.java}

about title: "GATK exome pipeline."
def BASEROOTDIR="BASEDIRNAME"

~~~~

Our pipeline makes several assumptions: 1) the protocol used is Illumina's TruSeq Exome Protocol, 2) the data is human, and 3) the raw files are named in the format: XXX_XXX_BARCODE_LANE_READ.fastq.gz.

~~~~{.java}

def PLATFORM="illumina"

~~~~

Resources files that are used throughout the GATK pipeline are declared as variables here so that can be reused throughout the pipeline; these files are from the GATK data bundle.

~~~~{.java}

def GENOME="$BASEROOTDIR" + "/bundle/2.8/hg19/ucsc.hg19.fasta"
def MILLS1000GINDEL="$BASEROOTDIR" +"/bundle/2.8/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
def DBDNP138="$BASEROOTDIR" +"/bundle/2.8/hg19/dbsnp_138.hg19.vcf"
def TARGETREGIONS="$BASEROOTDIR" +"/misc/TruSeq-Exome-Targeted-Regions.bed"
def HAPMAP="$BASEROOTDIR" +"/bundle/2.8/hg19/hapmap_3.3.hg19.sites.vcf"
def OMNI="$BASEROOTDIR" +"/bundle/2.8/hg19/1000G_omni2.5.hg19.sites.vcf"

~~~~

The GATK pipeline uses several programs that are not part of the GATK. These programs are declared as variables here.

~~~~{.java}

def BWA ="$BASEROOTDIR" + "/programs/bwa-0.7.12/bwa"
def SAMTOOLS = "$BASEROOTDIR"+ "/programs/samtools-1.2/samtools"
def MARKDUPLICATES = "$BASEROOTDIR"+ "/programs/MarkDuplicates.jar"

~~~~

The exome protocol captures regions slightly outside of exons; we use an interval padding of 100 bp to reflect this and this parameter is embedded in the GATK variable.

~~~~{.java}

def GATK = "$BASEROOTDIR"+ "/programs/GenomeAnalysisTK.jar --interval_padding 100 "

~~~~

Two routines, namely ```get_id``` and ```get_sample```, are defined to obtain metadata directly from the raw FASTQ file names. The routines expect files to be in the format: XXX_XXX_BARCODE_LANE_READ.fastq.gz. These are then used to add read groups (RG) to the BAM files.

~~~~{.java}

def get_id(filename) {
   def info = filename.split("/")[-1].split("\\.")[0].split("_")
   def lane = info[-2]
   def id = info[-4].concat('.').concat(info[-2])
   return(id)
}

def get_sample(filename) {
   def info = filename.split("/")[-1].split("\\.")[0].split("_")
   return(info[0])
}

~~~~

A script, ```check_bam_files.sh```, is used to check that the correct read groups (RG) are added to each BAM file. Errors may be introduced during parallelisation, when variables that should be local to one thread are used globally in other threads, which leads to namespace corruption. Using ```def``` should prevent this but it is important to double check that correct RGs are added.

~~~~{.java}

def CHECKBAM = "$BASEROOTDIR"+ "/scripts/check_bam_files.sh"

~~~~

The first step of the pipeline is to align the FASTQ files to a reference; hg19 is used in this pipeline and BWA-MEM is used as the alignment program as recommended by the GATK best practices. This is the step where the read groups are added to the BAM file.

The transform() function is a convenient alias for produce(); it deduces the name of the output or outputs by using the name of the input(s) and by modifying the file extension according to the transform() input.

~~~~{.java}

bwa_mem_align = {
   transform("bam") {
      def ID = get_id(input)
      def SAMPLE = get_sample(input)
      doc "Align reads to reference using BWA mem"
      output.dir = "alignments"
      var seed_length : 19
      
      exec """
      $BWA mem -M -t $threads -k $seed_length
      -R "@RG\\tID:${ID}\\tPL:$PLATFORM\\tPU:None\\tLB:None\\tSM:${SAMPLE}"
      $GENOME $input1 $input2 |
      $SAMTOOLS view -F 0x100 -bSu - | $SAMTOOLS sort - ${output.prefix}
      ""","bwamem"
   }
}

~~~~

The routine below simply indexes a BAM file. A forward instruction overrides the default files that are used for inputs to the next pipeline stage with files that you specify explicitly. We still want to use the BAM files as input in the next step.

~~~~{.java}

index_bam = {
   doc "Create BAM file index"
   output.dir=file(input.bam).absoluteFile.parentFile.absolutePath
   transform("bam") to ("bam.bai") {
      exec "$SAMTOOLS index $input.bam","index_bam"
   }
   forward input
}

~~~~

This is the step that runs the script for checking whether read groups match the file name.

~~~~{.java}

check_bam_for_sample_names = {
   doc "Checks if read groups in a bam file match the file name."
   exec "$CHECKBAM $input.bam"
   forward input
}

~~~~

Duplicates are creating during the DNA preparation step (e.g. during PCR amplification) and can cause biases that skew variant calling results. We use Picard to [remove potential PCR duplicates](https://www.broadinstitute.org/gatk/events/slides/1409/GATKwr5-BP-1-Map_and_Dedup.pdf), which are easy to spot because duplicate reads have an identical starting position and CIGAR string.

~~~~{.java}

dedup = {
   doc "Remove PCR duplicate reads from BAM"
   output.dir="alignments"
   exec """
   java -Xmx6g -jar $MARKDUPLICATES
   INPUT=$input.bam
   REMOVE_DUPLICATES=true
   VALIDATION_STRINGENCY=LENIENT
   AS=true
   METRICS_FILE=$output.metrics
   OUTPUT=$output.bam
   """
}

~~~~

Insertions and deletions (INDEL) cause misalignments, leading to false positive single nucleotide variants. The [step below](https://www.broadinstitute.org/gatk/events/slides/1409/GATKwr5-BP-2-Realignment.pdf) identifies potential INDEL regions that require realignment.

~~~~{.java}

realignIntervals = {
   doc "Select regions for realignment with GATK in 'realign' stage"
   output.dir="alignments"
   exec """
   java -Xmx4g -jar $GATK
   -T RealignerTargetCreator
   -R $GENOME
   -L $TARGETREGIONS
   -I $input.bam
   -nt $threads
   --known $MILLS1000GINDEL
   -o $output.intervals
   """
}

~~~~

Once candidate regions have been identified in the realignIntervals step, these regions are carefully realigned.

~~~~{.java}

realign = {
   doc "Apply GATK local realignment to intervals selected in stage 'realignIntervals' "
   output.dir="alignments"
   exec """
   java -Xmx5g -jar $GATK
   -T IndelRealigner
   -R $GENOME
   -I $input.bam
   -targetIntervals
   $input.intervals
   -o $output.bam
   """
}

~~~~

The [Base Quality Score Recalibration](https://www.broadinstitute.org/gatk/events/slides/1409/GATKwr5-BP-3-Base_recalibration.pdf) (BQSR) step is used to adjust base quality scores assigned by the sequencer; this is necessary as these base qualities are inaccurate and biased. A sliding window approach is used to tally the number of mismatches against the reference; loci that are known to vary are not used in the calculation. BQSR identifies patterns in how errors correlate with base features, such as the machine cycle or positions in a read.

~~~~{.java}

recal_count = {
   doc "Recalibrate base qualities in a BAM file using observed error rates"
   output.dir="alignments"
   INDEL_QUALS=""
   exec """
   java -Xmx5g -jar $GATK
   -T BaseRecalibrator
   -I $input.bam
   -nct $threads
   -R $GENOME
   -L $TARGETREGIONS
   --knownSites $DBDNP138
   --knownSites $MILLS1000GINDEL
   -l INFO
   -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate
   -o $output.counts
   """
}

~~~~

This is the second part of the BQSR step, which applies the recalibrated data to a new file.

~~~~{.java}

recal = {
   doc "Apply recalibration quality adjustments"
   output.dir="alignments"
   exec """
   java -Xmx4g -jar $GATK
   -T PrintReads
   -nct $threads
   -I $input.bam
   -BQSR $input.counts
   -R $GENOME
   -l INFO
   -o $output.bam
   """
}

~~~~

This step merges BAM files, if they belong to the same sample. If there is only one BAM file, the BAM file is copied with a new file name. The produce statement declares a block of statements that will be executed transactionally to create a given set of outputs.

~~~~{.java}

bam_merge = {
   doc "Merge Bam files "
   output.dir="alignments"
   
   def len = inputs.size()
   def SAMPLE = get_sample(input)
   def output = "${SAMPLE}" + ".bam"
   produce (output){
      if(len == 1){
         exec """
         cp  $inputs $output1
         ""","bam_merge"
      } else {
         exec """
         $SAMTOOLS merge -@ $threads $output1 $inputs
         ""","bam_merge"
      }
   }
}

~~~~

The [HaplotypeCaller](https://www.broadinstitute.org/gatk/guide/article?id=4148) is composed of four steps: defining active regions, determining haplotypes by re-assembly of the active regions, determining the likelihoods of the haplotypes given the read data, and assigning sample genotypes. We use the HaplotypeCaller in the ```GVCF``` mode, which allows variants to be called individually on each sample but merged later.

~~~~{.java}

HaplotypeCaller = {
   doc "Runs GATK HaplotypeCaller"
   
   output.dir="variants"
   transform("bam") to ("raw.snps.indels.g.vcf") {
      exec """
      java -Xmx6g -jar $GATK
      -T HaplotypeCaller
      -R $GENOME
      -I $input
      -nct $threads
      --emitRefConfidence GVCF
      --variant_index_type LINEAR
      --variant_index_parameter 128000
      --dbsnp $DBDNP138
      --max_alternate_alleles 50
      -L $TARGETREGIONS
      -o $output
      """
   }
}

~~~~

This step simply creates a file containing the list of [GVCF](https://www.broadinstitute.org/gatk/guide/article?id=4017) files.

~~~~{.java}

makevariantlist = {
   output.dir="variants"
   doc "Makes a list of all variant files"
   
   produce("variant.list"){
      exec """
      find variants -name '*.raw.snps.indels.g.vcf' > $output
      ""","makevariantlist"
   }
}

~~~~

This step combines the list of GVCF files.

~~~~{.java}

combine_GVCF = {
   output.dir="variants"
   
   doc "Combines variant files"
   produce ("combined.g.vcf"){
      exec """
      java -Xmx6g -jar $GATK
      -R $GENOME
      -T CombineGVCFs
      -L $TARGETREGIONS
      --variant $input
      -o $output
      """
   }
}

~~~~

This step will produce genotype likelihoods and re-genotype the newly merged record from the combine_GVCF step, and then re-annotate it.

~~~~{.java}

GenotypeGVCFs = {
   output.dir="variants"

   doc "Perform joint genotyping on gVCF files produced by HaplotypeCaller"
   produce ("genotype_all.vcf"){
      exec """
      java -Xmx10g -jar $GATK
      -R $GENOME
      -T GenotypeGVCFs
      -L $TARGETREGIONS
      --max_alternate_alleles 50
      -nt $threads
      --variant $input
      -o $output
      """
   }
}

~~~~

This step is the first part of the Variant Quality Score Recalibration (VQSR) [process](https://www.broadinstitute.org/gatk/guide/article?id=2805), which is the process of assigning accurate confidence scores to each putative variant call. The [first pass](https://www.broadinstitute.org/gatk/events/slides/1409/GATKwr5-BP-5-Variant_recalibration.pdf) consists of creating a Gaussian mixture model by looking at the distribution of annotation values over a high quality subset of the input call set, and then scoring all input variants according to the model.

~~~~{.java}

recalibrate_SNP = {
   output.dir="variants"

   doc "Build a recalibration model to score variant quality for filtering purposes"
   transform("snp.recal.csv","snp.tranches","snp.plot.R"){
      exec """
      java -Xmx10g -jar $GATK -nt $threads  
      -T VariantRecalibrator
      -R $GENOME
      -input $input
      -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP
      -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI
      -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $DBDNP138
      -an QD
      -an FS
      -an MQ
      -an MQRankSum
      -an ReadPosRankSum
      -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0
      -recalFile $output1
      -tranchesFile  $output2
      -rscriptFile  $output3
      """
   }
   forward input.vcf
}

~~~~

This is the second part of the VQSR step and consists of filtering variants based on score cut-offs identified in the first pass. The filter command is a convenient alias for produce where the name of the output or outputs is deduced from the name of the input(s) by keeping the same file extension but adding a new tag to the file name.

~~~~{.java}

Apply_SNP_Recalibration = {
   output.dir="variants"
   filter("snp_recalibrated"){
      exec """
      java -Xmx10g -jar $GATK
      -nt $threads 
      -T ApplyRecalibration
      -R $GENOME
      -input $input.vcf
      -mode SNP
      --ts_filter_level 99.0 
      -recalFile $input.csv
      -tranchesFile $input.tranches
      -o $output
      """
   }
}

~~~~

This is step one of the VQSR step for Indels.

~~~~{.java}

recalibrate_INDEL = {
   output.dir="variants"
   transform("indel.recal.csv","indel.tranches","indel.plot.R"){
      exec """
      java -Xmx10g -jar $GATK -nt $threads 
      -T VariantRecalibrator 
      -R $GENOME 
      -input $input.vcf 
      -resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS1000GINDEL 
      -an QD
      -an FS
      -an MQ
      -an MQRankSum
      -an ReadPosRankSum
      -mode INDEL 
      -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0
      --maxGaussians 4
      -recalFile $output1
      -tranchesFile $output2 
      -rscriptFile $output3
      """
   }
   forward input.vcf
}

~~~~

This step two of the VQSR step for Indels.

~~~~{.java}

Apply_INDEL_Recalibration = {
   output.dir="variants"
   filter ("indel_recalibrated"){
      exec """
      java -Xmx10g -jar $GATK -nt 32 
      -T ApplyRecalibration 
      -R $GENOME   
      -input $input.vcf
      -mode INDEL 
      --ts_filter_level 99.0 
      -recalFile $input.csv
      -tranchesFile $input.tranches 
      -o  $output.vcf
      """
   }
}

~~~~

Finally the code below defines the pipeline and how it should be run. The pipeline starts with an input splitting pattern defined by the ```%``` character; this defines which part of the file name should be used to split the input files into groups. The ```*``` character orders the groups but DOES NOT split the input; for more information take a look at [parallelising tasks](https://code.google.com/p/bpipe/wiki/ParallelTasks) in the Bpipe documentation.

~~~~{.java}

Bpipe.run {"%_L00%_R*.fastq.gz" * [ bwa_mem_align + index_bam + dedup + index_bam + realignIntervals + realign + index_bam + recal_count + recal + index_bam] + "%_*.fastq.dedup.realign.recal.bam"  * [bam_merge + check_bam_for_sample_names + index_bam + dedup + index_bam + realignIntervals + realign + index_bam + HaplotypeCaller ] + makevariantlist + combine_GVCF + GenotypeGVCFs + recalibrate_SNP + Apply_SNP_Recalibration + recalibrate_INDEL + Apply_INDEL_Recalibration}

~~~~

# Running the pipeline

To run the pipeline:

~~~~{.bash}
bpipe run -n 20 -r gatk_bwa_illunmina.groovy *.fastq.gz
~~~~

The ```-r``` parameter generates a basic report and the ```-n``` parameter defines the number of threads to use throughout the pipeline.

If the pipeline completed successfully, the file ```genotype_all.snp_recalibrated.indel_recalibrated.vcf``` should reside in the ```variants``` folder.

