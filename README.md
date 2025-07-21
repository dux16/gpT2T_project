# gpT2T_project
**Custom scripts used for analyses of paper [ Huizhong Fan<sup>#</sup>, Chentao Yang<sup>#</sup>, Xin Du<sup>#</sup>, Lei Nie<sup>#</sup>, Jiazheng Jin<sup>#</sup>, Chen Huang, Xin Huang, Yang Zhou, Xupeng Bi, Mengyang Xu, Jiang Hu, Guiqiang Wang, Lingxin Qiu, Dongya Wu, Ye Wang, Rong Hou, Guojie Zhang<sup>\*</sup>, Fuwen Wei<sup>\*</sup> ](https://xxx)
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
  
  # additional scaffolding process using hic data with yahs
  genome='PolarBear.hifiasm_ul_hic.asm.hic.p_ctg.fasta'
  hic1='xxx/polar_bear/00.dataset/Fastq/hic/BJX-xin_Hic.1.clean.fastq.gz'
  hic2='xxx/polar_bear/00.dataset/Fastq/hic/BJX-xin_Hic.2.clean.fastq.gz'
  # mapping
  xxx/software/Aligner/bwa/bwa mem -t 24 ../$genome $hic1 $hic2 |samtools view -@4 -bS -h - > aligned.bam
  samtools sort -@28 -o aligned_sorted.bam aligned.bam
  samtools index aligned_sorted.bam
  
  # YAHS 
  xxx/software/Assembly/yahs/yahs --no-contig-ec --telo-motif "TTAGGG" -o PolarBear ../$genome aligned_sorted.bam
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
  a. PAR bundary && X-Degenerate
   #lastZ X-Y alignment
    lastZ chrX.fa chrY.fa --ungapped --filter=identity:80 --filter=nmatch:400 --hspthresh=36400 --format=general-:name1,start1,end1,name2,start2,end2,strand2,nmatch >GPchrX2chrY.anchors
    lastZ --segments=GPchrX2chrY.anchors --filter=identity:80 --filter=nmatch:1000 --allocate:traceback=800M --format=general:name1,zstart1,end1,name2,strand2,zstart2+,end2+,nmatch,length1,id%,blastid% --rdotplot+score=GPchrX2chrY.dots | perl -alne 'if($F[4] eq "-"){($F[1], $F[2]) = ($F[2], $F[1])} print join "\t", @F'  >GPchrX2chrY.txt
   #mummer X-Y alignment & SyRI structure variants
    nucmer chrX.fa chrY.fa
    delta-filter -m -i 90 -l 100 out.delta >out.filter.delta
    show-coords -THrd out.filter.delta >out.filter.coords
    syri -c out.filter.coords -d out.filter.delta -r chrX.fa -q chrY.fa
  b.Ampliconic
   #LastZ self-alignment
    lastz chrY_softmask.fasta chrY_softmask.fasta --ungapped --filter=identity:80 --filter=nmatch:400 --hspthresh=36400 --format=general-:name1,start1,end1,name2,start2,end2,strand2,nmatch > chrY2chrY.anchors
    lastz chrY_softmask.fasta chrY_softmask.fasta --segments=GPchrY2chrY.anchors --filter=identity:80 --filter=nmatch:1000 --allocate:traceback=800M --format=general:name1,zstart1,end1,name2,strand2,zstart2+,end2+,nmatch,length1,id%,blastid% --rdotplot+score=chrY2chrY.dots > chrY2chrY.dat
    perl -lane 'if($F[4] eq "-"){($F[1], $F[2]) = ($F[2], $F[1])} $out=join("\t", @F); print $out' chrY2chrY.dat > chrY2chrY.dat.mod
    Rscript dotplot.X_ref.r chrY_softmask.fasta.fai chrY_softmask.fasta.fai chrY2chrY.dat.mod chrY_class.bed chrY2self_lastz_dotplot.pdf
  #blastn  self-alignment
    bedtools maskfasta -fi chrY.fa -bed chrY_hardmask.bed -fo chrY.hardmask.fasta
    seqtk subseq chrY.hardmask.fasta $MSY > MSY_hardmask.fasta
    makeblastdb -in MSY_hardmask.fasta -parse_seqids -dbtype nucl -out db/MSY
    samtools faidx MSY_hardmask.fasta
    awk -v OFS='\t' {'print $1,$2'} MSY_hardmask.fasta.fai > genome.txt
    bedtools makewindows -g genome.txt -w 5000 -s 2000 -i srcwinnum > genome.windows.w5ks2k.bed
    bedtools getfasta -fi MSY_hardmask.fasta -bed genome.windows.w5ks2k.bed > genome.windows.w5ks2k.fasta
    blastn -query genome.windows.w5ks1k.fasta -db db/MSY -outfmt 6 -perc_identity 50 -num_threads 32 -out merged_w5ks2k.blastn.outfmt6.txt
    perl -alne 'print if($F[3]>=2000)' merged_w5ks1k.blastn.outfmt6.txt > merged_w5ks2k.blastn.outfmt6.flt.txt
    Rscript createBedFileFromBlastn.R merged_w5ks2k.blastn.outfmt6.flt.txt
    bedtools merge -i merged_w5ks2k.blastn.outfmt6.flt.bed -d 1 -c 4 -o mean > merged_w5ks2k.blastn.outfmt6.flt.merged.bed
    bedGraphToBigWig merged_w5ks2k.blastn.outfmt6.flt.merged.bed chrY.len merged_w5ks2k.blastn.outfmt6.flt.merged.bw
  ```
#### - Composite Element Identification
  ```
  a. Moddotplot:
   for i in {5,10,20,40,60,80,100,120,140,150};do Wind=` exprt $i * 1000 `; moddotplot static -c config.json -w $Wind -o modeplot_${i}k;done
  b. Composite element
    ###annotation Composite element
     cat gpT2T_v1.0.genome.gff3 gpT2T_chrY_v1.0_final_pseudogenes.gff3  |awk '$3=="mRNA" && $1=="chrY"' >chrY_gene_pseodogene.gff
     awk '$1=="chrY"' GiantPanda.filteredRepeats.gff >chrY_repeat.gff
     cat chrY_gene_pseodogene.gff chrY_repeat.gff | awk '($4>=16262172 && $5<=20643784)||($4>=23249031 && $5<=39199034)||($4>=40012893 && $5<=45522213)||($4>=46436979)' |grep -E -w "LINE/L1|Satellite|Unknown|LTR/ERV1|mRNA"|sort -k4,4n >SDR_repeat_gene.gff
     awk '$3=="mRNA"' SDR_repeat_gene.gff |awk '{$1"\t"$4-1"\t"$5}' |sort -k2,2n |bedtools merge -i - -d 30000|awk '{print $1"\t"$2-50000"\t"$3+50000}'  >merge.region.bed
     awk '{if($3!="mRNA")print $1"\t"$4-1"\t"$5"\t"$7"\t"$3}' SDR_repeat_gene.gff >tmp1.bed
     awk -F "[\t=;]" '{if($3=="mRNA")print $1"\t"$4-1"\t"$5"\t"$7"\t"$NF}' >tmp2.bed
     cat tmp1.bed tmp2.bed |sort -k2,2n >tmp3.bed
     bedtools intersect -a tmp3.bed -b merge.region.bed -wa -wb >merge.region.intersect.bed
     python Class_Composite_element_mod.py merge.region.intersect.bed model.stat.bed
    ####Composite element phylogeny
     awk -F "[\t:-] "'{print $1"\t"$2"\t"$3}' model.stat.bed >region.bed
     bedtools getfasta -fi chrY.fa -bed region.bed >region.bed.fa
     samtools faidx chrY.fa chrY:470063453-47127767 >ref.fa
     ViralMSA.py -s region.bed.fa -r ref.fa -o viral_output -e test
     trimal -sgc -in region.bed.removeRef.fasta.aln >region.bed.removeRef.fasta.aln.stat
     trimal -gt 0.5 -in region.bed.removeRef.fasta.aln> region.bed.removeRef.fasta.trimal_0.5.aln.fasta
     raxml-ng -msa region.bed.removeRef.fasta.aln.rename.trimal_0.5.fasta --model GTR+G
  ```
#### - ERV1 Analysis
  ```
  n=0
  grep ERV1 gpT2T_earlGrey_Rep_annot.bed |grep -v chrY |awk '{if($3-$2>=2000)print $1"\t"$2+1"\t"$3}'|while read a b c;do n=`expr $n + 1`;echo -e "$a\t$b\t$c\t${a}_erv$n";done >ERV1.region.txt
  grep ERV1 gpT2T_earlGrey_Rep_annot.bed |grep  chrY >chrY.erv.bed
  cat chrY.erv.bed | awk '{if(($2>=16262177 && $3<=20643784)||($2>=23249031 && $3<=39199034) ||($2>=40012893 && $3>=45522213)|| ($2>=46435704 && $3<=47471033))print}' >chrY.L1_erv_sat.bed
  n=0
  cat chrY.erv.bed chrY.L1_erv_sat.bed |sort |uniq -c |awk '{if($4-$3>=2000)print}' |while read a b c d e f g;do n=`expr $n + 1`;c=`expr $c + 1` if [ $a -gt 1 ];then echo -e "$b\t$c\t$d\t${b}_erv$n" ;else echo -e "$b\t$c\t$d\t${b}L1_erv$n";fi;done >>ERV1.region.txt
  cat ERV1.region.txt |while read a b c d ;do seq=`samtools faidx gpT2T_v1.0.genome.fna $a:$b-$c|grep -v ">"`;echo -e ">$d\n$seq";done >ERV1.fa
  samtools faidx ERV1.fa chr1_erv12  >reference.fa
  mafft --auto --thread 10 ERV1.fa >ERV1.fa.align
  raxml-ng --msa ERV1.fa.align --model GTR+G
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
