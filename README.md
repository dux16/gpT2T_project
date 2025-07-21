# gpT2T_project
**Custom scripts used for analyses of paper [ Huizhong Fan<sup>#</sup>, Chentao Yang<sup>#</sup>, Xin Du<sup>#</sup>, Lei Nie<sup>#</sup>, Jiazheng Jin<sup>#</sup>, Chen Huang, Xin Huang, Yang Zhou, Xupeng Bi, Mengyang Xu, Jiang Hu, Guiqiang Wang, Linxin Qiu, Dongya Wu, Ye Wang, Rong Hou, Guojie Zhang<sup>\*</sup>, Fuwen Wei<sup>\*</sup> ](https://xxx)
.**<br>

In this study, we did not develop new software; thus, we provide example commands and in-house scripts for analyses or provide links to other sources employed.

### Content:
 - [Genome Assembly](https://github.com/dux16/gpT2T_project#--genome-assembly)<br>
   * [Initial assembly](https://github.com/dux16/gpT2T_project#--contig-assembly)<br>
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
   * [Centromere Prediction](https://github.com/dux16/gpT2T_project#--centromere-prediction)<br>
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
  yak count -k31 -b37 -t10 -o father.yak panda_father_80x_R1.fastq.gz panda_father_80x_R2.fastq.gz
  yak count -k31 -b37 -t10 -o mother.yak panda_mother_80x_R1.fastq.gz panda_mother_80x_R2.fastq.gz
  ### Perfoming hifiasm
  hifiasm -1 father.yak -2 mother.yak -o gpT2T.hifiasm -t 80 $path/m64*.ccs.fq.gz
  gfatools gfa2fa gpT2T.hifiasm.dip.hap1.p_ctg.gfa > gpT2T.hifiasm.dip.hap1.p_ctg.fna
  gfatools gfa2fa gpT2T.hifiasm.dip.hap2.p_ctg.gfa > gpT2T.hifiasm.dip.hap2.p_ctg.fna
  
  ## Verkko assemblies combining HiFi, Ultra-Long ONT, and NGS reads
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
  ### HicTrioBinning
  bash ./HTB.sh -M gpT2T_mat.fasta -P gpT2T_pat.fasta -1 GP_HiC_R1.fastq.gz -2 GP_HiC_R2.fastq.gz -N 100 -B bwa -S seqtk -O ./ -I 1
  cat maternal.reads_1.fq.gz homo.reads_1.fq.gz > mat_HiC_R1.fastq.gz
  cat maternal.reads_2.fq.gz homo.reads_2.fq.gz > mat_HiC_R2.fastq.gz
  cat paternal.reads_1.fq.gz homo.reads_1.fq.gz > pat_HiC_R1.fastq.gz
  cat paternal.reads_2.fq.gz homo.reads_2.fq.gz > pat_HiC_R2.fastq.gz
  ### Juicer and  3D-DNA pipeline for maternal
  bwa index assembly.haplotype1.fasta
  python2 $juicer/misc/generate_site_positions.py DpnII mat assembly.haplotype1.fasta
  
  awk 'BEGIN{OFS="\t"}{print $1, $NF}' mat_DpnII.txt > mat.chrom.sizes
  bash ./scripts/juicer.sh -g gpT2T_mat -z assembly.haplotype1.fasta -p mat.chrom.sizes -y mat_DpnII.txt -D /public/home/duxin/software/05.3D_genome/juicer-1.6 -d /data01/gpasm/01.GP_T2T_asm/04.hap_Chr_asm/02.chr_asm/mat -t 60 -S early
  ### Juicer and  3D-DNA pipeline for paternal
  bwa index assembly.haplotype2.fasta
  python2 $juicer/misc/generate_site_positions.py DpnII pat assembly.haplotype2.fasta
  
  awk 'BEGIN{OFS="\t"}{print $1, $NF}' pat_DpnII.txt > pat.chrom.sizes
  bash ./scripts/juicer.sh -g gpT2T_pat -z assembly.haplotype2.fasta -p pat.chrom.sizes -y pat_DpnII.txt -D /public/home/duxin/software/05.3D_genome/juicer-1.6 -d /data01/gpasm/01.GP_T2T_asm/04.hap_Chr_asm/02.chr_asm/pat -t 60 -S early
  ```
  The [JuiceBox](https://s3.us-east-1.wasabisys.com/hicfiles/public/Juicebox/Juicebox_1.11.08.exe) was used to manually curate the Hi-C maps for maternal and paternal.

#### - Gap filling
  ```
  ## Making the haplotype-specific 21-mer markers
  meryl k=21 count output mat_R1.meryl mat_80x_R1.fastq.gz
  meryl k=21 count output mat_R2.meryl mat_80x_R2.fastq.gz
  meryl k=21 count output pat_R1.meryl pat_80x_R1.fastq.gz
  meryl k=21 count output pat_R2.meryl pat_80x_R2.fastq.gz
  meryl k=21 count output child_R1.meryl child_R1.fq.gz
  meryl k=21 count output child_R2.meryl child_R2.fq.gz
  meryl union-sum output mat.meryl mat_R*.meryl
  meryl union-sum output pat.meryl pat_R*.meryl
  meryl union-sum output child.meryl child_R*.meryl
  $MERQURY/trio/hapmers.sh mat.meryl/ pat.meryl child.meryl
  meryl print mat_only.meryl | awk '{print $1}' > maternal.uniq.flt.mer
  meryl print pat_only.meryl | awk '{print $1}' > paternal.uniq.flt.mer
  ## Running HAST4TGS: https://github.com/BGI-Qingdao/HAST4TGS
  CLASSIFY_ONLY.sh --paternal_mer paternal.uniq.flt.mer --maternal_mer maternal.uniq.flt.mer --format fastq --thread 80 --offspring 20220406-UNL272-P4-PAK05241.pass.fastq.gz --offspring 20220406-UNL272-P4-PAK05745.pass.fastq.gz --offspring 20220406-UNL272-P4-PAK05782.pass.fastq.gz --offspring 20220406-UNL272-P4-PAK05820.pass.fastq.gz --offspring 20220406-UNL272-P4-PAK05887.pass.fastq.gz --offspring 20220407-UNL272-P4-PAK06109.pass.fastq.gz --offspring 20220407-UNL272-P4-PAK06244.pass.fastq.gz --offspring 20220407-UNL272-P6-PAK05314.pass.fastq.gz --offspring 20220407-UNL272-P6-PAK05756.pass.fastq.gz --offspring 20220407-UNL272-P6-PAK05829.pass.fastq.gz --offspring 20220407-UNL272-P6-PAK07349.pass.fastq.gz --offspring 20220407-UNL272-P6-PAK07797.pass.fastq.gz
  ## Applying new standard
  perl renew_UL-ONT_trio-binning.pl ./phasing.out
  ## Extracting UL-ONT reads
  for fq in `ls *pass.fq.gz`; do for nf in maternal_new.cut paternal_new.cut; do hap=${nf/_new.cut/}; p=${fq/.fastq.gz/}; seqtk subseq $fq $nf | pigz -9 --best -p 20 > ${hap}_$p.fq.gz; done; done
  ## minimap2
  
  ## GapCloser: https://github.com/BGI-Qingdao/TGS-GapCloser
  ```
#### - rDNA assembly
  ```
  ## Running rnammer to identify 18S and 28S rRNA
  rnammer -S euk -gff mat_rRNA.gff -xml mat_rRNA.xml gpT2T.hifiasm.dip.hap1.p_ctg.fna
  rnammer -S euk -gff pat_rRNA.gff -xml pat_rRNA.xml gpT2T.hifiasm.dip.hap2.p_ctg.fna
  ## Running infernal to identify 5.8S rRNA
  cmscan -Z length*2/1000000 --cut_ga --rfam --nohmmonly --tblout mat_tempout/seq1.tblout --fmt 2 --cpu 20 --clanin /xxx/infernal-1.1.4-linux-intel-gcc/Rfam_lib/Rfam.clanin /xxx/infernal-1.1.4-linux-intel-gcc/Rfam_lib/Rfam.cm mat-split/seq1.fa > mat_tempout/seq1.cmscan
  cmscan -Z length*2/1000000 --cut_ga --rfam --nohmmonly --tblout pat_tempout/seq2.tblout --fmt 2 --cpu 20 --clanin /xxx/infernal-1.1.4-linux-intel-gcc/Rfam_lib/Rfam.clanin /xxx/infernal-1.1.4-linux-intel-gcc/Rfam_lib/Rfam.cm pat-split/seq2.fa > pat_tempout/seq2.cmscan
  perl convert_cmscan2gff3.pl mat_tempout mat_infernal
  perl convert_cmscan2gff3.pl pat_tempout pat_infernal
  ## Mummer
  mummer -maxmatch -l 20 -b -F -L -c rDNA_region.fasta rDNA_region.fasta > rDNA_region2self.out
  mummerplot -p rDNA_region2self rDNA_region2self.out --png
  ## Identifying rDNA units of maternal and paternal assemblies
  
  ## Detecting the copy number of rDNA in maternal and paternal genomes based on Trio-Bining UL-ONT reads
  
  ```
#### - Genome polishing and assessment
  ```
  ## PacBio HiFi polishing
  
  ```
#### - Polar bear genome assembly
  ```
  date
  # Hybrid assembly with HiFi, ultralong and Hi-C reads
  hic1=xxx/polar_bear/00.dataset/Fastq/hic/BJX-xin_Hic.1.clean.fastq.gz
  hic2=xxx/polar_bear/00.dataset/Fastq/hic/BJX-xin_Hic.2.clean.fastq.gz
  ul=xxx/polar_bear/00.dataset/Fastq/nanopore/ONT_50k.fastq.gz,xxx/polar_bear/00.dataset/Fastq/nanopore/added.ONT_50k.fastq.gz
  hifi="xxx/polar_bear/00.dataset/Fastq/hifi/HiFi.fastq.gz xxx/polar_bear/00.dataset/Fastq/hifi/added.HiFi.fastq.gz"
  xxx/software/hifiasm-0.24.0/hifiasm -t 24 -o PolarBear.hifiasm_ul_hic.asm --h1 $hic1 --h2 $hic2 --ul $ul $hifi
  date
  
  ```

### - Genome Annoation
#### - Annotation of repetitive elements
  ```
  ## 01. RepeatMasker
  ## 02. RepeatModeler
  ## 03. RepeatProteinMasker
  ## 04. TRF
  ```
#### - Annoation of protein coding genes
  ```
  ## Homology-based prediction
  ## RNA-seq annotation
  ## Ab initio prediction
  ## EVM
  ```

### - Genomic Analyses
#### - Heterozygous Variant Analyses
  ```
  ## for SNP, indel detection
  
   for i in `seq 1 20`
   do
   	chr="chr"$i
  	 [ -d $chr ] || mkdir $chr
  	 cd $chr
  	 echo "#!/bin/bash
   date
   $soft/nucmer --maxmatch -t 48 -l 100 -c 500 -p $chr $data/$chr/$chr.mat.fasta $data/$chr/$chr.pat.fasta
   $soft/delta-filter -m -i 90 -l 100 $chr.delta > $chr.delta.filt
   $soft/dnadiff -d $chr.delta.filt -p $chr
   mv $chr.snps $chr.var
   python3 scripts/split_mummerVar2snpAndindel.py $chr.var $chr
   python3 scripts/indel_statistic.py $chr.indel
   python3 scripts/small_indel_Pos2bed.py  $chr.indel.report.txt > $chr.indel.bed
   date" >  $chr.sh
  	 easy_submit  -t 48 -m 20G -s yes -p bnode $chr.sh
  	 cd ..
   done
   
  ## combine all chrs
  # snp
  for chr in `ls |grep chr|grep -v 'chrX' |grep -v 'chrY'`
  do
  	cat $chr/$chr.snp | awk '{print $11"\t"$1-1"\t"$1"\t"$2"\t"$3"\t"$12"\t"$4-1"\t"$4}'
  done > mat_vs_pat.snp.bed
  
  # small indel
  for chr in `ls |grep chr`
  do
  	cat $chr/$chr.indel.bed
  done  > mat_vs_pat.indel.bed
  
  ## for SV detection
  mumdir='xxx/nucmer_MatRef' # which contains the nucmer mapping results
  for a in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20
  do
  	[ -d $a  ] || mkdir $a
  	cd $a
  	[ -L $a.delta.filt  ] || ln -s $mumdir/$a/${a}.delta.filt
  	mat=`head -1 $a.delta.filt|awk '{print $1}'`
  	pat=`head -1 $a.delta.filt|awk '{print $2}'`
  
  	echo "#file	name	tags
  $mat	${a}_hap1	lw:1.5
  $pat	${a}_hap2	lw:1.5" > genome.txt
  	echo "#!/bin/bash
  date
  export PATH=\"xxx/software/Aligner/mummer4/bin:\$PATH\"
  show-coords -THrd ${a}.delta.filt > ${a}.coords.filt
  source xxx/miniconda3/bin/activate syri_env
  syri -c ${a}.coords.filt -d ${a}.delta.filt -r $mat -q $pat
  plotsr --sr syri.out --genomes genome.txt -H 1 -W 12 -o $a.plotsr.pdf
  echo '$a done'
  date" > ${a}.syri.sh
  
      sh ${a}.syri.sh
  	plotsr --sr syri.out --genomes genome.txt -H 1 -W 12 -o $a.plotsr.pdf
  	cd ..
  done
  
  find ./ -name "syri.out" | grep -v 'sex' | xargs cat |sort -k1,1V -k2n > autosome.syri.out
  python3 ../../../bin/post_syri.py autosome.syri.out > autosome.syri.sv.tab
  cut -f 12 autosome.syri.sv.tab|awk '{sum+=$1};END{print sum}'
  ```
#### - SNP Calling, Nucleotide Diversity and ROH Analyses
  ```
  
  ```
#### - Estimation of *De Novo* Mutations
  ```
  
  ```

### - Centromere Analysis
#### - Centromere prediction
  ```
  
  ```
#### - HOR Detection
  ```
  ```
#### - Composit Elements
  ```
  ```
### - Y Chromosome Analysis
#### - Y chromosome structure analysis
  ```
  
  ```
#### - Composite Element Identification
  ```
  ```
#### - ERV1 Analysis
  ```
  ```
#### - Curation of Y chromosome genes
   ```
   ```
#### - The species-species divergence calculation of giant panda and brown bear
  ```
  ```
#### - The analysis of autosomal transposition
  ```
  ```
#### - Gene conversion detection
  ```
  ```
#### - The transcriptome analysis of Y chromosome genes
  ```
  ```
### - Comparative Genomic Analysis
#### - Orthologous Genes Identification
  ```
  ```
#### - Selection Pressure Analysis
  ```
  ```
#### - Conserved Non-coding Elements Analysis
  ```
  ```
#### - Pseudogene Identification
  ```
  ```
#### - Functional Enrichment Analysis
  ```
  ```
