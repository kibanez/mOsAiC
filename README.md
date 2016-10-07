# mOsAiC

Somatic mosaicism brings up to the presence at the same time, in a single individual, of two genetically distinct cell populations arising from a postzygotic mutation.

Numerous computational approaches have been developed in order to detect somatic mosaic variants from NGS data. They usually do the analysis from paired tumor-normal samples, as Strelka (Saunders et al. 2012) or MuTect (Cibulskis et al. 2013). Strelka uses a Bayesian strategy representing continuous allele frequencies for both tumor and normal samples, and influencing the expected genotype of the normal. It assumes the normal samples as a mixture of germline variation with noise, and the tumor samples as a mixture of the normal sample with somatic variation (Saunders et al. 2012). The approach presented in MuTect also uses a Bayesian classifier in the variant detection in the tumor sample, once the low quality sequences have been removed, filters in order to remove false positives resulting from correlated sequencing artifacts not captured by the error model, and appoints to the variants as somatic or germline by a second Bayesian classifier (Cibulskis et al. 2013). These two strategies have described benchmarking approaches that use real, rather than simulated, sequencing data to evaluate the accuracy of each tool. However, due to the homology of several genomic locus within the panel design, low quality reads that could be mapped in other homologous regions are removed with these both strategies.


We here present an in-house developed approach in order to detect variants in mosaicism from a list of already aligned reads files (BAM files). The summary of the NGS preprocessing is as follows:

The patient samples were merged in the same pool, and sequenced. This detail is important to the post-processing step, since working with several samples is advantageous. For instance, it permits to identify variants that are recurrent along the samples in the pool, due to frequent polymorphism or sequencing artifacts.
Pair-end raw reads coming from Illumina NextSeq500 were converted into paired FASTQ files using bcl2fastq-v2.15.0.4 software from Illumina (https://github.com/brwnj/bcl2fastq). Then preprocessing was done by using Trimmomatic (Bolger et al. 2014), trimming and cropping FASTQ data as well as removing adapters.. Afterwards, balanced reads were mapped to the published human genome reference (UCSC hg19) by using Bowtie2 aligner (Langmead and Salzberg 2012). Then PCR duplicate reads were removed by using one of the most commonly used approach Picard MarkDuplicates (http://broadinstitute.github.io/picard/). Following the best practices proposed by GATK (McKenna et al. 2010; DePristo et al. 2011; Van der Auwera et al. 2013; Vrijenhoek et al. 2015), an INDEL realignment was done afterwards in order to determine suspicious intervals which were likely in need of realignment (GATK RealignerTargetCreator function) and then, a local realignment of reads was done to correct misalignments due to presence of INDELs (GATK IndelRealigner function). The input of this mosaicism tool we here develope are the aligned, realigned, and subsequent recalibrated reads (BAM files).



Step 1. mpileup
The first step is to extract the base-pair information for each genomic position from the BAM files, facilitating the subsequent SNP/INDEL calling. For such aim, samtools mpileup v1.3 (Li et al. 2009; Li 2011a) was used for each BAM file along the sequenced pool. The base quality (q value) and the mapping quality scores (Q value) were set to 0, since many genes in the capture design have pseudo-genes or high homology with other genomics sites. Thus, using recommended scores (q = 20, Q = 20) would include out reads mapped in these genes with low mapping quality scores defined by the aligner considering that they also map in other genomic position.

Step 2. Call capture
Once all the base-pair information was available for each sample, the variant capture calling was proceed. For such purpose, bcftools v1.3 (Li 2011) was applied. It is important to take into account, that both samtools and bcftools must be from the same version. This calling was done keeping all possible alternate alleles at each variant site (option -A) as well as keeping all the alternative models for multiallelic and rare-variant callings (option -m) in a VCF file. As consequence, in this step a VCF file was obtained containing the capture of the variant calling, considering the previous two conditions.

Step 3. AVAF computation
We here define a new attribute encompassing the alternative variant allele fraction (AVAF). Our objective is to detect and categorize all the alternative alleles in each genomic position along the design panel, with a minimum value of ratio of 0.02 (2%). We assume that alternative alleles with a frequency less than 0.02 might be due to background noise, error during the variant calling process or attributable to a low complex region.
The idea is to compute the ratio (i.e., percentage) in which each alternative variant allele appears in each genomic position, and label it depending on whether the AVAF is higher or lower than the mosaicism threshold given (m option in the tool). Particularly, the threshold used in this study was setup to 0.07. 
The purpose of this mosaicism tool was to empower the final user when analyzing the resulting entire table. Thus, this approach does not filter out any variant in which AVAF is less than m. Accordingly, when the AVAF is lower than m, the variant is labelled as LowDepthAltMosaicism. In contrast, whether AVAF is higher than m, the variant is identified as PASS. This classification states in the FILTER field in the corresponding output VCF file.

Step 4. Split multi-allelic sites into simple sites
In the VCF file there are genomic positions with an unique alternative allele, and other position in which multiple alternative alleles are detected. For a proper analysis, it is convenient and necessary to separate multi allelic sites into individual sites.
After this step, the resulting VCF file could contain many repetitive genomic positions with same reference alleles but different alternative alleles information. In this manner it is possible to study separately each alternative allele with its corresponding AVAF, among other attribute values.

Step 5. Variant recurrence determination
Once the AVAF attribute is defined for each alternative allele along all the genomic positions defined, the recurrence of each variant throughout all the samples analyzed with the same designed genomic regions was computed. This value was included in a new attribute called samplesRun. 
This step is essential since it permits to observe in how many samples examined with the same genomic bed file is detected each variant. Very high value of samplesRun might indicate very polymorphic variants or sequencing bias among others. In consequence, low samplesRun values (i.e., samplesRun <= 2) are very interesting, since it suggests that only in the studying sample or in another one is detected that variant.
In this step the percentage or ratio of the variant in a set of pseudo-control samples (11 healthy exomes) computation is also done.

Step 6. Germinal vs. somatic variants comparison
The contrasting of matched tumor (i.e., tissue) and germline (i.e., blood) samples is crucial to distinguish somatic from germline variants, bringing out low allelic fraction from background noise caused by the high-sequencing error rate. Hence, working with single or individual tumor samples would lead to an identification of a larger number of false positive events.
This step is independent from the previous part. Here the tissue and the blood VCF files from the same sample are studied, and also the VCF files from two healthy samples as control. If there is no a corresponding blood sample, only the tissue and the two control samples are examined. The process starts with all the variants detected in the tissue VCF. All the variants identified as PASS in the blood VCF are filtered out, since they are also detected in the tissue VCF and cannot be somatic mosaic variants. In contrast, variants identified as LowDepthAltMosaicism are not filtered out, but enriched with a new attribute called lowMosaic, which contains the corresponding AVAF value of the variant in the blood sample. When processing the control samples the procedure is stringent. Variants labelled as PASS or LowDepthAltMosaicism in the control samples that are also detected in the tissue VCF are separated out.


Two configuration (*cfg) example files are included: 

> configuration_mosaic_computation_part1.cfg configuration file compresses the Step 1 to Step5: the mosaic variant computation.
> configuration_somatic_mosaic_computation_part2.cfg configuration file compresses the Step 6. This step requires the previus step : the somatic mosaic variant computation.

NOTE: this tool does not annotate the resulting VCF files. This must be done independently.
