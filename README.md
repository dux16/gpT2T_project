# gpT2T_project
**Custom scripts used for analyses of paper [ Fan et al. xxx ](https://xxx)
.**<br>

In this study, we did not develop new software; thus, we provide example commands for analyses or provide links to other sources employed.

### Content:
 - [Genome Assembly](https://github.com/dux16/gpT2T_project#--genome-assembly)<br>
   * [Initial assembly)](https://github.com/dux16/gpT2T_project#--contig-assembly)<br>
   * [Pseudochromosome construction](https://github.com/dux16/gpT2T_project#--contig-assembly)<br>
   * [Gap filling](https://github.com/dux16/gpT2T_project#--polishing)<br>
   * [rDNA assembly](https://github.com/dux16/gpT2T_project#--rDNA-assembly)<br>
   * [Genome polishing and Assessment](https://github.com/dux16/gpT2T_project#--genome-polishing-and-assessment)<br>
   * [Polar bear genome assembly](https://github.com/dux16/gpT2T_project#--polar-bear-genome-assembly)<br>
   * [gpT2T vs GPv1](https://github.com/dux16/gpT2T_project#--gpT2T-vs-GPv1)<br>
 - [Genome Annotation](https://github.com/dux16/gpT2T_project#--genome-annotation)<br>
   * [Annoation of repetitive elements](https://github.com/dux16/gpT2T_project#--annotation-of-repetitive-elements)<br>
   * [Prediction of protein coding genes](https://github.com/dux16/gpT2T_project#--prediction-of-protein-coding-genes)<br>
 - [Genomic Analysis](https://github.com/dux16/gpT2T_project#--genomics-analysis)<br>
   * [Heterozygous Variant Analyses](https://https://github.com/dux16/gpT2T_project#--genome-annotation)<br>
   * [SNP Calling, Nucleotide Diversity and ROH Analyses](https://https://github.com/dux16/gpT2T_project#--SNP-calling-nucleotide-diversity-and-ROH-analyses)<br>
   * [Estimation of *De Novo* Mutations](https://github.com/dux16/gpT2T_project#--haplotype-analysis)<br>
 - [Centromere Analysis](https://github.com/dux16/gpT2T_project#--centromere-analysis)<br>
   * [Centromere Prediction]()<br>
   * [HOR Detection](https://github.com/dux16/gpT2T_project#--HOR-detection)<br>
   * [Composit Elements](https://github.com/dux16/gpT2T_project#--composit-elements)<br>
 - [Y Chromosome Analysis](https://github.com/dux16/gpT2T_project#--Y-chromosome-analysis)<br>
   * [Y chromosome structure analysis](https://github.com/dux16/gpT2T_project#--Y-chromosome-structure-analysis)<br>
   * [Composite Element Identification](https://github.com/dux16/gpT2T_project#--composite-element-identification)<br>
   * [ERV1 Analysis](https://github.com/dux16/gpT2T_project#--ERV1-analysis)<br>
   * [Curation of Y chromosome genes](https://github.com/dux16/gpT2T_project#--curation-of-Y-chromosome-genes)<br>
   * [The species-species divergence calculation of giant panda and brown bear](https://github.com/dux16/gpT2T_project#--the-species-species-divergence-calculation-of-giant-panda-and-brown-bear)<br>
   * [The analysis of autosomal transposition](https://github.com/dux16/gpT2T_project#--the-analysis-of-autosomal-transposition)<br>
   * [Gene conversion detection](https://github.com/dux16/gpT2T_project#--gene-conversion-detection)<br>
   * [The transcriptome analysis of Y chromosome genes](https://github.com/dux16/gpT2T_project#--the-transcriptome-analysis-of-Y-chromosome-genes)<br>
 - [Comparative Genomic Analysis](https://github.com/dux16/gpT2T_project#--comparative-genomics-analysis)<br>
   * [Orthologous Genes Identification](https://github.com/dux16/gpT2T_project#--orthologous-genes-identification)<br>
   * [Selection Pressure Analysis](https://github.com/dux16/gpT2T_project#--selection-pressure-analysis)<br>
   * [Conserved Non-coding Elements Analysis](https://github.com/dux16/gpT2T_project#--conserved-non-coding-elements-analysis)<br>
   * [Pseudogene Identification](https://github.com/dux16/gpT2T_project#--pseudogene-identification)<br>
   * [Functional Enrichment Analysis](https://github.com/dux16/gpT2T_project#--functional-enrichment)<br>


## Analyses
### - Genome assembly

#### - Initial assembly
  ```
  ## HiFiasm assemblies
  ### Building the 31-mer hash table of panda parents using yak
  yak count -k31 -b37 -t10 -o father.yak panda_father_80x_1.fastq panda_father_80x_2.fastq
  yak count -k31 -b37 -t10 -o mother.yak panda_mother_80x_1.fastq panda_mother_80x_2.fastq
  ### Perfoming hifiasm
  hifiasm -1 father.yak -2 mother.yak -o gpT2T.hifiasm -t 80 ./m64*.ccs.fq.gz
  gfatools 
  
  ## Verkko assemblies combining HiFi and Ultra-Long ONT reads
  ### Preparing the compressed homopolymer reads
  dehomopolymerate -f mat_80x_merged.fq.gz | pigz --fast -p 12 -9 > mat_nohomop.fa.gz
  dehomopolymerate -f pat_80x_merged.fq.gz | pigz --fast -p 12 -9 > pat_nohomop.fa.gz
  dehomopolymerate -f child_reads.fq.gz | pigz --fast -p 12 -9 > child_nohomop.fa.gz
  ### Running meryl and hapmers.sh to build 51-mer homopolymer-compressed database
  meryl k=51 count output mat_compress.meryl mat_nohomop.fa.gz
  meryl k=51 count output pat_compress.meryl pat_nohomop.fa.gz
  meryl k=51 count output child_compress.meryl child_nohomop.fa.gz
  sh $MERQURY/trio/hapmers.sh mat_compress.meryl pat_compress.meryl child_compress.meryl
  ### Running verkko pipeline
  verkko -d gpT2T --hifi $path/02.HiFi_reads/*gz --nano $path/01.ONT_reads/*gz --threads 40 --hap-kmers mat.hapmer.meryl pat.hapmer.meryl trio --sto-run 1 50 24 --mer-run 40 200 200 --ovb-run 40 150 200 --red-run 40 200 200 --mbg-run 40 200 200 --ali-run 40 200 200 --ruk-run 40 200 200 --sub-run 40 200 200 --cns-run 40 200 200 --par-run 40 200 200 --pbs --min-ont-length 10000
  
  ```
#### - Pseudochromosome construction
  ```
  ## Scaffolding
  ### [HicTrioBinning](https://github.com/BGI-Qingdao/HicTrioBinning)<br>
  bash ./HTB.sh -M gpT2T_mat.fasta -P gpT2T_pat.fasta -1 GP_HiC_R1.fastq.gz -2 GP_HiC_R2.fastq.gz -N 100 -B bwa -S seqtk -O ./ -I 1
  cat maternal.reads_1.fq.gz homo.reads_1.fq.gz > mat_HiC_R1.fastq.gz
  cat maternal.reads_2.fq.gz homo.reads_2.fq.gz > mat_HiC_R2.fastq.gz
  cat paternal.reads_1.fq.gz homo.reads_1.fq.gz > pat_HiC_R1.fastq.gz
  cat paternal.reads_2.fq.gz homo.reads_2.fq.gz > pat_HiC_R2.fastq.gz
  ### Juicer and  3D-DNA pipeline for maternal
  bwa-mem2 index assembly.haplotype1.fasta
  python2 $juicer/misc/generate_site_positions.py DpnII mat assembly.haplotype1.fasta
  
  awk 'BEGIN{OFS="\t"}{print $1, $NF}' mat_DpnII.txt > mat.chrom.sizes
  bash ./scripts/juicer.sh -g gpT2T_mat -z assembly.haplotype1.fasta -p mat.chrom.sizes -y mat_DpnII.txt -D /public/home/duxin/software/05.3D_genome/juicer-1.6 -d /data01/gpasm/01.GP_T2T_asm/04.hap_Chr_asm/02.chr_asm/mat -t 60 -S early
  ### Juicer and  3D-DNA pipeline for paternal
  ```
#### - Gap filling
  ```
  
  ```
