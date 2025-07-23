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
 - [Chromosome Evolution Analysis](https://github.com/dux16/gpT2T_project#--chromosome-evolution-analysis)<br>
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
  yak count -k31 -b37 -t10 -o father.yak pat_NGS_80x_R1.fastq.gz pat_NGS_80x_R2.fastq.gz
  yak count -k31 -b37 -t10 -o mother.yak mat_NGS_80x_R1.fastq.gz mat_NGS_80x_R2.fastq.gz
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
  meryl k=21 count output mat_R1.meryl mat_NGS_80x_R1.fastq.gz
  meryl k=21 count output mat_R2.meryl mat_NGS_80x_R2.fastq.gz
  meryl k=21 count output pat_R1.meryl pat_NGS_80x_R1.fastq.gz
  meryl k=21 count output pat_R2.meryl pat_NGS_80x_R2.fastq.gz
  meryl k=21 count output child_R1.meryl child_NGS_80x_R1.fq.gz
  meryl k=21 count output child_R2.meryl child_NGS_80x_R2.fq.gz
  meryl union-sum output mat.meryl mat_R*.meryl
  meryl union-sum output pat.meryl pat_R*.meryl
  meryl union-sum output child.meryl child_R*.meryl
  $MERQURY/trio/hapmers.sh mat.meryl/ pat.meryl child.meryl
  meryl print mat_only.meryl | awk '{print $1}' > maternal.uniq.flt.mer
  meryl print pat_only.meryl | awk '{print $1}' > paternal.uniq.flt.mer
  ## Running HAST4TGS: https://github.com/BGI-Qingdao/HAST4TGS
  CLASSIFY_ONLY.sh --paternal_mer paternal.uniq.flt.mer --maternal_mer maternal.uniq.flt.mer --format fastq --thread 80 --offspring 20220406-UNL272-P4-PAK05241.pass.fastq.gz --offspring 20220406-UNL272-P4-PAK05745.pass.fastq.gz --offspring 20220406-UNL272-P4-PAK05782.pass.fastq.gz --offspring 20220406-UNL272-P4-PAK05820.pass.fastq.gz --offspring 20220406-UNL272-P4-PAK05887.pass.fastq.gz --offspring 20220407-UNL272-P4-PAK06109.pass.fastq.gz --offspring 20220407-UNL272-P4-PAK06244.pass.fastq.gz --offspring 20220407-UNL272-P6-PAK05314.pass.fastq.gz --offspring 20220407-UNL272-P6-PAK05756.pass.fastq.gz --offspring 20220407-UNL272-P6-PAK05829.pass.fastq.gz --offspring 20220407-UNL272-P6-PAK07349.pass.fastq.gz --offspring 20220407-UNL272-P6-PAK07797.pass.fastq.gz
  ## Applying new standard
  perl ./scripts/renew_UL-ONT_trio-binning.pl ./phasing.out
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
  perl ./scripts/convert_cmscan2gff3.pl mat_tempout mat_infernal
  perl ./scripts/convert_cmscan2gff3.pl pat_tempout pat_infernal
  ## Mummer
  mummer -maxmatch -l 20 -b -F -L -c rDNA_region.fasta rDNA_region.fasta > rDNA_region2self.out
  mummerplot -p rDNA_region2self rDNA_region2self.out --png
  ## Identifying rDNA units of maternal and paternal assemblies
  
  ## Detecting the copy number of rDNA in maternal and paternal genomes based on Trio-Bining UL-ONT reads
  
  ```
#### - Genome polishing and assessment
  ```
  ## PacBio HiFi polishing
  nextpolish2 xxx
  ## Evaluating the quality of assemblies
  ### Building a merged 21-mer database
  meryl k=21 count output child_NGS1.meryl child_NGS_R1.fq.gz
  meryl k=21 count output child_NGS1.meryl child_NGS_R1.fq.gz
  meryl k=21 count output child_hifi.meryl gp_HiFi.fastq.gz
  meryl union-sum output merged.meryl child_*.meryl
  ### Running Merqury
  /xxx/merqury/merqury.sh merged.meryl gpT2T_mat_polished.genome.fna gpT2T_pat_polished.genome.fna gpT2T_polish
  ### Manually select the chromosome with the highest quality value to generate the chimeric genome
  chimeric=gpT2T.genome.fna 
  /xxx/merqury/merqury.sh merged.meryl $chimeric gpT2T
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
  #### - gpT2T vs GPv1
  ```
  
  ```
### - Genome Annoation
#### - Annotation of repetitive elements
  ```
  ## 01. RepeatMasker
  RepeatMasker -pa 50 -nolow -norna -no_is -a -gff -species mammals gpT2T.genome.fna
  perl ./scripts/ConvertRepeatMasker2gff.pl gpT2T.genome.fna.out TE.gff TE
  ## 02. RepeatModeler
  BuildDatabase -name gpT2T gpT2T.genome.fna
  RepeatModeler -pa 50 -database gpT2T -LTRStruct
  RepeatMasker -pa 50 -lib ./custom.lib gpT2T.genome.fna -a
  perl ./scripts/ConvertRepeatMasker2gff.pl gpT2T.genome.fna.out Denovo.gff Denovo
  ## 03. RepeatProteinMasker
  perl ./scripts/split_fasta.pl gpT2T.genome.fna split_by_scaffold 20000000
  for i in split_by_scaffold/*.fa; do echo "RepeatProteinMask -noLowSimple -pvalue 1e-04 -engine rmblast $i; perl ./scripts/mod_RPM_out.pl $i.annot > $i.annot.mod";  done > 02.split.running.sh
  ParaFly -CPU 50 -c 02.split.running.sh
  cat split_by_scaffold/*annot.mod > gpT2T.genome.fna.repeatproteinmasker.annot
  rm -rf split_by_scaffold
  perl ./scripts/ConvertRepeatMasker2gff.pl gpT2T.genome.fna.repeatproteinmasker.annot TP.gff TP
  ## 04. TRF
  perl ./scripts/split_fasta.pl gpT2T.genome.fna split_by_scaffold 5000000
  for i in split_by_scaffold/*.fa; do i=`basename $i`; echo "cd `pwd`/split_by_scaffold; timeout 7200 trf $i 2 7 7 80 10 50 2000 -d -h; perl mod_trf_out.pl ${i}.2.7.7.80.10.50.2000.dat > ${i}.2.7.7.80.10.50.2000.dat.mod"; done > 02.split_run_trf.sh
  ParaFly -CPU 50 -c 02.split_run_trf.sh
  find ./split_by_scaffold -name "*.dat.mod" | xargs cat  > gpT2T.genome.fna.2.7.7.80.10.50.2000.dat
  rm -rf split_by_scaffold
  perl ./scripts/ConvertTrf2Gff.pl gpT2T.genome.fna.2.7.7.80.10.50.2000.dat TRF.gff
  ```
#### - Annoation of protein coding genes
  ```
  ## Homology-based prediction
  for i in cat cattle dog giant_panda horse human pig polar_bear;
  do
  echo "java -Xmx150g -jar /xxx/GeMoMa-1.9.jar CLI GeMoMaPipeline threads=50 AnnotationFinalizer.r=NO pc=true o=true t=/xxx/chimeric_genome/gpT2T.genome.fna GAF.a=\"iAA>=0.75\" i=$i a=/xxx/z.homology_species/${i}.genome.gff g=/xxx/z.homology_species/${i}.genome.fna outdir=/xxx/$i GeMoMa.m=500000; perl ./scripts/ConvertFormat_GeMoMa.pl /xxx/$i/final_annotation.gff"
  done > run_GeMoMa.sh
  echo "
#SBATCH --job-name GMM
#SBATCH -c 50 --mem 150G
#SBATCH --get-user-env
#SBATCH -o /xxx/GMM.slurm.out
#SBATCH -e /xxx/GMM.slurm.err
#SBATCH -D /xxx
#SBATCH -N 1
#SBATCH --partition=bnode
bash run_GeMoMa.sh
" > GeMoMa_sbatch.sh
  sbatch GeMoMa_sbatch.sh
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
  ## Running pav to call heterozygous variants
  singularity run pav_latest.sif  -c 48
  ## Regions of high heterozygosity
  echo "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tVARIANTS/KB" > hetero_50kb.density
  bedtools coverage -a <(bedtools makewindows -g gpT2T_v1.0.genome.fa.fai) -w 50000 | awk '{OFS="\t";print $1,$2+1,$3}') -b pav_gpT2T.vcf.gz -counts | awk 'BEGIN{scale=10;OFS="\t"}{print $1,$2,$3,$4,$4/50}'>> hetero_50kb.density
  ## Regions of high nucleotide diversity
  sort -k5gr gp.50kb.pi.windowed.pi|head -n 4862 | sort -k1V | cut -f1,4,5,9 -d$'\t' > gp.50kb_0.1top.bed
  bedtools intersect -a gpT2T_annot_v2.0_PCGs.gff3 -b gp.50kb_0.1top.bed -u -wa | cut -f1,4,5,9 -d$'\t' | awk '{for (i=1; i<=NF; i++) if ($i ~ /Name=/) {sub(/;.*$/, "", $i); split($i, a, "="); print a[2]}}' > gp.50kb_0.1top_genes.list
  ## ROHs
  multiIntersectBed -i ./00.data/QIN*.ROH.bed |awk '{OFS="\t";print $1,$2,$3,$3-$2,$4}' > QIN_isecSTAT_temp.bed
  multiIntersectBed -i ./00.data/SC*.ROH.bed |awk '{OFS="\t";print $1,$2,$3,$3-$2,$4}' > SC_isecSTAT_temp.bed
  multiIntersectBed -i ./00.data/*.ROH.bed |awk '{OFS="\t";print $1,$2,$3,$3-$2,$4}' > ISEC_isecSTAT_temp.bed
  awk '$5>=8 {OFS="\t";print $1,$2,$3}' QIN_isecSTAT_temp.bed > QIN_ROH_union_indiv8_temp.bed
  awk '$5>=17 {OFS="\t";print $1,$2,$3}' SC_isecSTAT_temp.bed > SC_ROH_union_indiv17_temp.bed
  bedtools intersect -a QIN_ROH_union_indiv8_temp.bed -b SC_ROH_union_indiv17_temp.bed > QIN_SC_ROH_union_isec_indiv_temp.bed
  bedtools subtract -a QIN_ROH_union_indiv8_temp.bed -b SC_ROH_union_indiv17_temp.bed > QIN_ROH_union_sub_indiv_temp.bed
  bedtools subtract -b QIN_ROH_union_indiv8_temp.bed -a SC_ROH_union_indiv17_temp.bed > SC_ROH_union_sub_indiv_temp.bed
  echo -e "CHROM\tSTART\tEND\tindivCNT\tPROPORTION\tLABEL" > result.bed
  bedtools intersect -a <(cut -f1-3,5 ISEC_isecSTAT_temp.bed) -b QIN_SC_ROH_union_isec_indiv_temp.bed | uniq | awk 'BEGIN{OFS="\t"}{print $0,$4/48,"ISEC"}' >> result.bed
  bedtools intersect -a <(cut -f1-3,5 QIN_isecSTAT_temp.bed) -b QIN_ROH_union_sub_indiv_temp.bed | uniq | awk 'BEGIN{OFS="\t"}{print $0,$4/15,"QIN"}' >> result.bed
  bedtools intersect -a <(cut -f1-3,5 SC_isecSTAT_temp.bed) -b SC_ROH_union_sub_indiv_temp.bed | uniq | awk 'BEGIN{OFS="\t"}{print $0,$4/33,"SC"}' >> result.bed
  sort -k1,1V -k2,2n result.bed -o result.bed && rm *_temp*
  bedtools intersect -a gpT2T_annot_v2.0_PCGs.gff3 -b <(grep QIN result.bed) -u -wa | awk '{OFS="\t";print $1,$4,$5,$9}' | awk '{for (i=1; i<=NF; i++) if ($i ~ /Name=/) {sub(/;.*$/, "", $i); split($i, a, "="); print a[2]}}' > QIN_sub_genes.list
  ```
#### - Estimation of *De Novo* Mutations
  ```
  
  ```

### - Centromere Analysis
#### - Centromere prediction
  ```
  mkdir tmp

## run trf to find 
trf gp_assembly.fasta 4 7 7 80 10 500 2000 -d -h


## transform dat to bed
awk 'BEGIN{chr=""; OFS="\t"} {if($1 == "Sequence:") {chr=$2} else if(NF==15){print chr,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}}' gp_assembly.fasta.4.7.7.80.10.500.2000.dat | awk 'OFS="\t"{print $1,$2,$3,"repeat_"NR,$9,$4,$5,$6,$7,$8,$10,$11,$12,$13,$14,$15}'> pandas.over5k.bed 

#filter region with length <= 5k or similarity < 70 or motif length < 20 or containing less than 20 motifs and generate motif fasta
bedmap --max-element --fraction-either 0.1 pandas.over5k.bed | uniq > pandas.over5k.no_overlap.bed
awk '$7>20 && $6>= 20 && $9 >= 70 && $3-$2 > 5000' pandas.over5k.no_overlap.bed > pandas.over5k.no_overlap.idover70.bed
sort pandas.over5k.no_overlap.idover70.bed | uniq | awk '{print ">"$4"\n"$16}' > pandas.over5k.no_overlap.idover70.pri.fasta

## remove similar motifs
cd-hit-est -i pandas.over5k.no_overlap.idover70.pri.fasta -o pandas.over5k.no_overlap.idover70.fasta -d 0 -aS 0.95 -c 0.95 -G 1 -g 1

sort pandas.over5k.no_overlap.idover70.bed | uniq | awk '{print ">"$4"\n"$16""$16}' > pandas.over5k.no_overlap.idover70.double.fasta
minimap2 -xmap-ont -c -k 11 -p 0.2 pandas.over5k.no_overlap.idover70.double.fasta pandas.over5k.no_overlap.idover70.fasta | awk '$1 != $6 && $11/$2>0.8' | sort -k2,2n -k1,1 -k6,6> self_align.paf
cut -f1,6 self_align.paf | uniq > pair.ids
python3 Partition.py -rel pair.ids | cut -f2- | sed 's/\t/\n/g' > dup.ids 
python3 fastaKit.py -exl dup.ids pandas.over5k.no_overlap.idover70.fasta > pandas.over5k.no_overlap.idover70.rmdup.fasta 
python3 fastaKit.py -exl dup.ids pandas.over5k.no_overlap.idover70.double.fasta > pandas.over5k.no_overlap.idover70.rmdup.double.fasta
minimap2 -xmap-ont -c -k 11 -p 0.2 pandas.over5k.no_overlap.idover70.rmdup.double.fasta pandas.over5k.no_overlap.idover70.rmdup.fasta | awk '$1 != $6 && $11/$2>0.9' | sort -k2,2n -k1,1 -k6,6> self_align.r2.paf
cut -f1,6 self_align.r2.paf | uniq > pair.r2.ids 
python3 Partition.py -rel pair.r2.ids | cut -f2- | sed 's/\t/\n/g' > dup.r2.ids
python3 fastaKit.py -exl dup.r2.ids pandas.over5k.no_overlap.idover70.rmdup.fasta > pandas.over5k.no_overlap.idover70.rmdup.r2.fasta 
python3 fastaKit.py pandas.over5k.no_overlap.idover70.rmdup.r2.fasta -len > pandas.over5k.no_overlap.idover70.rmdup.r2.len
trf pandas.over5k.no_overlap.idover70.rmdup.r2.fasta  2 7 7 80 10 30 2000 -d -h > repeat.fasta
mv pandas.over5k.no_overlap.idover70.rmdup.r2.fasta.2.7.7.80.10.30.2000.dat tmp/pandas_repeat.trf.dat
awk 'BEGIN{chr=""; OFS="\t"} {if($1 == "Sequence:") {chr=$2} else if(NF==15){print chr,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}}' tmp/pandas_repeat.trf.dat | awk 'OFS="\t"{print $1,$2,$3,"subrepeat_"NR,$9,$4,$5,$6,$7,$8,$10,$11,$12,$13,$14,$15}'> tmp/pandas_repeat.trf.bed
cut -f1 tmp/pandas_repeat.trf.bed | sort | uniq > ./tmp/exist.ids
grep -wf ./tmp/exist.ids pandas.over5k.no_overlap.idover70.rmdup.r2.len > tmp/pandas_repeat.trf.len 
join tmp/pandas_repeat.trf.len tmp/pandas_repeat.trf.bed | awk '($4-$3)/$2 > 0.9' | awk 'OFS="\t" {print $1,$3,$4,$5,$6/($4-$3),$17}' > tmp/pandas_repeat.trf_contain_sub.bed
bedmap --max-element --fraction-either 0.1 tmp/pandas_repeat.trf_contain_sub.bed | sort |uniq > tmp/pandas_repeat.trf_contain_sub.rm_overlap.bed
cut -f1 tmp/pandas_repeat.trf_contain_sub.rm_overlap.bed > tmp/pandas_repeat.trf_contain_sub.rm_overlap.ids
python3 fastaKit.py -exl tmp/pandas_repeat.trf_contain_sub.rm_overlap.ids pandas.over5k.no_overlap.idover70.rmdup.r2.fasta > tmp/pandas_repeat.rmdup.r2.pure.fasta
sort tmp/pandas_repeat.trf_contain_sub.rm_overlap.bed | uniq | awk '{print ">"$1"\n"$6}' > tmp/pandas_repeat.trf_contain_sub.rm_overlap.fasta
python3 fastaKit.py -lgt 20 tmp/pandas_repeat.trf_contain_sub.rm_overlap.fasta > tmp/pandas_repeat.trf_contain_sub.rm_overlap.long.fasta
cat tmp/pandas_repeat.rmdup.r2.pure.fasta tmp/pandas_repeat.trf_contain_sub.rm_overlap.long.fasta > tmp/pandas_repeat.rmdup.r2.rmsub.fasta
awk '{if($1 ~ ">") {print $0} else {print $0$0}}' tmp/pandas_repeat.rmdup.r2.rmsub.fasta > tmp/pandas_repeat.rmdup.r2.rmsub.double.fasta 
minimap2 -xmap-ont -c -k 11 -p 0.2 tmp/pandas_repeat.rmdup.r2.rmsub.double.fasta tmp/pandas_repeat.rmdup.r2.rmsub.fasta | awk '$1 != $6 && $11/$2>0.9' | sort -k2,2n -k1,1 -k6,6> tmp/self_align.r3.paf
cut -f1,6 tmp/self_align.r3.paf | uniq > tmp/pair.r3.ids 
python3 Partition.py -rel tmp/pair.r3.ids | cut -f2- | sed 's/\t/\n/g' > tmp/dup.r3.ids
python3 fastaKit.py -exl tmp/dup.r3.ids tmp/pandas_repeat.rmdup.r2.rmsub.fasta > pandas.over5k.no_overlap.idover70.rmdup.rmsub.r3.fasta

## re-align the motif to assembly
RepeatMasker -pa 25 -cutoff 100 -lib pandas.over5k.no_overlap.idover70.rmdup.rmsub.r3.fasta -s -nolow -gff -dir repeatmasker gp_assembly.fasta
  ```
#### - HOR Detection
  ```
## generate hicat HOR annotation. hor_region.bed contain records of the Sat_556 arrays and mask.bed contain records of non-Sat_556 motifs (include TE and composite element) annotations, mon.fasta is a SAT_556 motif
cat "hor_region.bed" | while read a; do
    chr=`echo $a | cut -f 1 -d " "`
    begin=`echo $a | cut -f 2 -d " "`
    end=`echo $a | cut -f 3 -d " "`
    dirname="${chr}"
    [ -d $dirname ] || mkdir $dirname
    cd $dirname
    mkdir tmp
echo "#!/bin/bash
grep -w ${chr} ../mask.bed > mask.bed
python3 getInterval.py -all -mask --cutfile mask.bed -fasta ${fastadic}/${chr}.fasta > masked.fa
python3 getInterval.py -c -low ${begin} -up ${end} -fasta masked.fa | sed 's/_/-/g'> ${chr}.HOR.fasta
python3 HiCAT.py -i ${chr}.HOR.fasta -t ../../mon.fasta -o ${chr}_hor -st 0.01 -th 18
" > ${chr}.work.sh
    sbatch -c 18 --mem 30000 --partition=cpu ${chr}.work.sh
    cd ..
done

## generate NTRprism plots to revalidate the HOR annotation. The NTRprism script can be found at https://github.com/altemose/NTRprism
python3 getInterval.py -all -mask --cutfile mask.bed -fasta ~/pandas_T2T/data/gp_mat.v0.4.6.polish2.fasta > masked.fa
python3 cutbyBed.py -fasta masked.fa -bed hor_region.bed > hor_region.fasta
perl NTRprism_ProcessFasta_v0.22.pl hor_region.fasta test 1 100000 30 17 0
for a in `ls test.*.span100000.k17.mincount30.bin1.txt`; do 
filename=`basename $a`
name=`echo $filename | sed -e 's/.span100000.k17.mincount30.bin1.txt//g' -e 's/test.region_//g'`
chrom=`echo $name | cut -f1 -d "."`
pos=`echo $name | cut -f2,3 -d "."`
grep -v NNNNNNNNNNNNNNNNN $a > ${chrom}.span100000.k17.mincount30.bin1.txt
Rscript NTRprism_PlotSpectrum.r --args chrX.span100000.k17.mincount30.bin1.txt $pos 10000 $chrom
done
  ```
#### - Composit Elements
  ```
  ## dotplot generate via gepard
java -cp Gepard-2.1.jar org.gepard.client.cmdline.CommandLine -seq ${chr}.fasta -matrix edna.mat -maxwidth 1250 -maxheight 1250 -word 50 -from1 ${start} -to1 ${end} -from2 ${start} -to2 ${end} -format bmp -outfile ${chr}_50.bmp

  ## The unannotated regions are found via dotplot and the following region is used as motif to run Repeatmasker
  ## chr1_hap2	106808338	106890009
  ## chr3_hap2	87393465	87670972
  ## chr17_hap2	13358594	14460267
  ## chr18_hap2	43963293	44006075
  ## chr20_hap2	27092662	27525483
  RepeatMasker -pa 10 -lib gp_assembly.fasta -s -gff -dir repeatlib ./motif.potential.fasta

  ## filter the resulting Repeatmasker out file. 
  awk 'OFS="\t"{if(NR>3)if(NF==16){print $0}else {print $0"."}}' gp_assembly.fasta.out | awk 'OFS="\t" {if($9=="+"){print $5,$6-1,$7,$10,100-$2,$9,$16} else {print $5,$6-1,$7,$10,100-$2,"-",$16}}' | sort -k4,4 | join -1 4 - motif.length | awk 'OFS="\t" {if($5 >= 95) print $2,$3,$4,$1,$5,$6,$7,$8,($4-$3)/$8}' | awk 'OFS="\t" {if($4 == "poten_4" && ($3-$2) > 7000){print $1,$2,$3,$4,$5,$6,$7,$8,$9} else if($4 == "poten_3" && $3-$2> 9000){print $1,$2,$3,$4,$5,$6,$7,$8,$9} else if($9>0.4 && $4 != "poten_4" && $4 != "poten_3") {print $1,$2,$3,$4,$5,$6,$7,$8,$9}}'> potential_comp.len_over80.bed
  awk 'OFS="\t" {print $1":"$4,$2,$3}' potential_comp.len_over80.bed | sort -k1,1 -k2,2n |bedtools merge -i - -d 1500 | awk 'OFS="\t"{print $0,$3-$2}' | sed 's/:/\t/g' | awk 'OFS="\t" {print $1,$3,$4,$2,$4-$3}'> test.merge.bed
  awk 'OFS="\t"{if(NR>3)if(NF==16){print $0}else {print $0"."}}' gp_mat.v0.4.6.polish2.fasta.out | awk 'OFS="\t" {if($9=="+"){print $5,$6-1,$7,$10,100-$2,$9,$16,$12,$13,$14} else {print $5,$6-1,$7,$10,100-$2,"-",$16,$14,$13,$12}}' | sort -k4,4 -k1,1V > potential_comp.all.bed
  bedtools subtract -a potential_comp.all.bed -b gp_mat.rp_r1-r2.bed -f 0.8 -A | awk '$5 > 90 '| sort -k4,4 -k1,1V -k2,2n > check.bed

  ## filter comp_TE1 Repeatmasker result seperately.
  awk 'OFS="\t"{if(NR>3)if(NF==16){print $0}else {print $0"."}}' gp_assembly.fasta.out | awk 'OFS="\t" {if($9=="+"){print $5,$6-1,$7,$10,100-$2,$9,$12,$13,$16} else {print $5,$6-1,$7,$10,100-$2,"-",$14,$13,$16}}' | sort -k4,4 | join -1 4 - motif.length | awk 'OFS="\t" {if($5>90 && $1 == "poten_4") print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10}' > potent4_extend.bed
  sort -k1,1V -k2,2n potent4_extend.bed | bedtools merge -i - -d 2000 | awk '$3-$2> 10000' |bedtools subtract -a - -b potential_comp.region.back.bed -N -f 0.9 | awk 'OFS="\t" {print $1,$2,$3,"poten4_"NR,$3-$2}' | bedtools intersect -a - -b ../../../r2/centromere_region/mat.cen_ext.v1.bed -wa > poten4.addition.bed 
  bedtools intersect -a poten4.addition.bed -b potent4_extend.bed -wo | awk 'OFS="\t"{print $4,$12,$13}' > motif.info.bed
  cut -f1 motif.info.bed | sort | uniq |awk 'OFS="\t" {print $1,"0","10742"}' > motif.siminfo.bed
  bedtools coverage -a motif.siminfo.bed -b motif.info.bed | awk '$7 > 0.8' | cut -f1 > conserve.ids
  grep -wf conserve.ids poten4.addition.bed > poten4.check.bed

  ## Final interval is manully determined from check.bed and poten4.check.bed
  ```
### - Chromosome Evolution Analysis
#### Infer the ancestral karyotype using DESCHRAMBLER

> #### Input file preparation

```bash
# In the early stage, it is necessary to run Lastz to obtain high-quality Chain and Net files;
# The chromosome number and the corresponding file name can only contain lowercase "chr";
mkdir DESCHRAMBLER && cd DESCHRAMBLER
mkdir chainNet
#The directory structure of chainNet is the reference species AME, which contains directories of other species. The subdirectories are chain and net, for example:
chainNet/
└── AME   #reference species
    ├── CLU   # query species
    │   ├── chain
    │   │   ├── chr1.chain  # <Lastz_path>/4.chaincleaner/all_clean.chain
    │   │   ├── ...
    │   │   └── chrX.chain
    │   └── net
    │       ├── chr1.net   # <Lastz_path>/7.netfilter/target.filtered.net
    │       ├── ...
    │       └── chrX.net
    ├── FCA   # query species
#Copy chain and Net files
mkdir -p chainNet/AME/CLU/chain chainNet/AME/CLU/net
for i in `seq 1 20` X; do cp /work1/fanhzh/project/gp_project/2.karyotype/1.lastz/Chr$i/*CLU.simple.genome/4.chaincleaner/all_clean.chain chainNet/AME/CLU/chain/chr${i}.chain; cp /work1/fanhzh/project/gp_project/2.karyotype/1.lastz/Chr$i/*CLU.simple.genome/7.netfilter/target.filtered.net chainNet/AME/CLU/net/chr${i}.net; done

#tree.txt (timetree,The function of @ is to mark ancestral branches)
vi tree.txt
(((SSU:36.96320000,(FCA:11.91318000,PTI:11.91318000):25.05002000):18.39844000,((VLA:12.24414000,CLU:12.24414000):32.85586000,(((MFL:13.50000000,(LLU:10.83500000,(MER:6.80000000,NVI:6.80000000):4.03500000):2.66500000):26.62000000,(ZCA:23.80995000,NSC:23.80995000):16.31005000):0.10000000,(AME:19.25625000,(TOR:13.8381800,(UAR:0.72000000,UMA:0.72000000):13.1181800):5.31807000):20.86375000):4.98000000):10.26164000)@:18.88836000,MPE:74.25000000);

#Prepare three files: config.SFs, Makefile.SFs and params.txt, as required (https://github.com/jkimlab/DESCHRAMBLER)
```

> #### Run DESCHRAMBLER

```bash
perl ~/Software/DESCHRAMBLER/DESCHRAMBLER.pl ./params.txt >DESCHRAMBLER.log 2>&1
```

> #### Visualization

```bash
python -m jcvi.graphics.chromosome AME.feature.bed APCF.idmap --size=AME.sizes
```

```bash
# cat AME.feature.bed
# Deal with APCF_AME.merged.map, reformatted as follows:
chr18   13233350        19663033        APCF1
chr18   53277   12616899        APCF1
chr18   19768604        25532631        APCF1
```

```bash
# cat AME.sort.sizes
chr10   114930382
chr11   116659357
chr1    203573159
```

```bash
# cat APCF.idmap
APCF1   1       #C56673
APCF2   2       #D66B61
APCF3   3       #DD7D58
APCF4   4       #E3AB51
APCF5   5       #EACC89
APCF6   6       #E5DE59
APCF7   7       #D2D553
APCF8   8       #A1C363
APCF9   9       #62B16B
APCF10  10      #45A16F
APCF11  11      #488B6C
APCF12  12      #6BB795
APCF13  13      #45ABA6
APCF14  14      #4A8AAC
APCF15  15      #545F99
APCF16  16      #805697
APCF17  17      #A94B7B
APCF18  18      #CD5689
APCF19  19      #D1CBB2
APCF20  20      #BF9E78
APCF21  21      #AD8A6E
```
> #### The way to obtain xxxx.feature.bed and xxxx.sizes
```bash
grep APCF APCF_AME.merged.map > 1; grep chr APCF_AME.merged.map > 2; paste 2 1 > 3; awk '{print $1"\t"$2"\t"$3}' 3 > 4; sed -i "s/[.:-]/ /g" 4; sed -i "s/\t \t/\t-\t/g" 4 ; awk '{print $2"\t"$3"\t"$4"\t"$6$7"\t"$5}' 4 > AME.feature.bed; rm 1 2 3 4;
awk '{print $1"\t"$2}' ../../../01.rawdata/genome/AME.genome.fna.fai > AME.sizes; sed -i "s/Chr/chr/g" AME.sizes
```
> #### Reconstruct the ancestral chromosome sequence based on the correspondence of the bed file
```perl
perl Construct_Ancestor_Genome.pl AME.feature.bed AME.genome.fna Ancestor.fasta
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
> ##### Whole genome alignment
# Whole genome alignment pipeline

> ##### 1. Before running the process, it is necessary to simplify the processing of genomic files and GFF files, extract chromosomes and simplify chromosome names(Chr*), and extract the longest transcript.

> ##### 2. Preparation of the directory and scripts before the whole genome alignment.

```bash
faSplit byname ./gpT2T_v1.0.genome.fna . #split ref genome
for i in `seq 1 20` X Y; do mkdir Chr$i; mv Chr$i.fa Chr$i; done #Adjust the parameters here according to the chromosome quantity

for i in `seq 1 20` X Y; do cp lastz_Chr.sh Chr$i/; sed -i "s/Chr1/Chr$i/g" Chr$i/lastz_Chr.sh; done # Batch rename chromosome IDs in lastz_Chr.sh under each folder

for i in `seq 1 20` X Y; do echo "cd Chr$i && sh lastz_Chr.sh && sh command_Chr$i.list && cd .."; done >command_lastz.sh

bash command_lastz.sh

for i in ./Chr*; do root=`basename $i`; cd $i && perl ../easy_submit.pl -ss SLURM -m 30G -t 10 -p long -s yes -l 4 command_${root}.list && cd ..; done   #To improve runtime efficiency, submit jobs via the SLURM system. If you don’t have a cluster, split the command_Chr*.list into multiple execution files by chromosome and run them in parallel using ParaFly

### sh this scrpits and modify the following script
```

> ##### 3. Multiz

```bash
#Running multiz
python runmutliz.py #revise path and parameter based on your data

for i in `seq 1 20` X Y; do echo "sh Chr$i/multiz.sh"; done >command_multiz.list

perl ../../easy_submit.pl -ss SLURM -m 15G -t 5 -p normal -s yes -l 1 command_multiz.list  #This step is for submitting commands to the cluster; if you don't have a cluster, you can skip it.

#ParaFly -c command_split_maf.list -CPU 22
```
> ##### Getting Orthologous genes
> ##### 2. MaftoGene

```bash
mkdir genes
mkdir align
for i in `seq 1 22` X Y
do
mkdir $i
echo "perl ./scripts/01.convertMaf2List.pl /xxx/1.CNE/1.multiz/Chr$i/tmp.18.maf $i/Chr$i.gene HSA AJU,AME,BTA,CLU,FCA,LCA,LLU,MER,MPE,NNE,NVI,PBE,PCR,PTI,SIN,UAM,UAR,UMA,VLA 4" >> maf2gene.sh
echo "perl ./scripts/02.lst2gene.pl $i/Chr$i.gene /xxx/HSA.genome.gff" >> lst2gene.sh
done

perl /xxx/easy_submit.pl -ss SLURM -m 15G -t 1 -p normal -s yes -l 1 maf2gene.sh

perl /xxx/easy_submit.pl -ss SLURM -m 15G -t 1 -p normal -s yes -l 1 lst2gene.sh
```

> ##### 3. Filtergenes

```bash
perl ./scripts/filter_genes_base_gap.pl genes gene HSA,BTA,MPE HSA 0.2
```

> ##### 4. MACSE alignment

```bash
ls genes/*.fas > genelist; sed -i 's/\.fas//g' genelist

mkdir align_AA align_NT

cat genelist | while read line; do echo "java -jar -Xmx20G /work1/fanhzh/project/gp_project/software/macse_v2.07.jar -prog alignSequences -seq align/${line}.fas -gap_op -7 -gap_ext -1 -fs -30 -gc_def 1 -stop -100 -out_AA align_AA/${line}_AA.fas -out_NT align_NT/${line}_NT.fas";done > run_alignment.sh

ParaFly -c run_alignment.sh -CPU 30
```

> ##### 5. MACSE export

```bash
mkdir stats align_AA_filter align_NT_filter

cat genelist | while read line; do echo "java -jar -Xmx600m /work1/fanhzh/project/gp_project/software/macse_v2.07.jar -prog exportAlignment -align align_NT/${line}_NT.fas -cons_threshold 0.7 -codonForFinalStop --- -codonForInternalStop NNN -codonForExternalFS --- -codonForInternalFS NNN -charForRemainingFS - -out_NT align_NT_filter/${line}_noFS_NT.fasta -out_AA align_AA_filter/${line}_noFS_AA.fasta -out_stat_per_seq stats/${line}_stats.csv"; done > exportAlignment.sh

ParaFly -c exportAlignment.sh -CPU 30
```

#### - Selection Pressure Analysis

> ##### 1. Reconstruct the phylogenetic tree.

```bash
mkdir con_tre
perl merge_singlecopygenes_phy.pl align_AA_filter con_tree
cd con_tree
perl ../scripts/paml_scripts/Fasta2Phylip.pl trim.merged.aln.fasta trim.merged.aln.phy

echo "raxmlHPC-PTHREADS-AVX2 -T 32 -f a -N 100 --bootstop-perms=1000 -m PROTGAMMAWAG -x 12345 -p 12345 -o "HSA,BTA,MPE" -s trim.merged.aln.phy -n trim.merged.aln.nwk" > contree.sh
perl ../../../easy_submit.pl -ss SLURM -m 50G -t 32 -p long -s yes contree.sh

#Prepare tree file based on the trim.merged.aln.nwk result file
cd ../ && mkdir paml && cd paml && mkdir tree && cd tree
vi tree.nwk
20 1
(HSA,BTA,(MPE,((((((FCA,PBE),AJU),LCA),(PTI,NNE)),PCR),((CLU,VLA),((((UAR,UMA),UAM),AME 1),(((NVI,MER),LLU),SIN))))));

vi tree.label.nwk
20 1
(HSA,BTA,(MPE,((((((FCA,PBE),AJU),LCA),(PTI,NNE)),PCR),((CLU,VLA),((((UAR,UMA),UAM),AME #1),(((NVI,MER),LLU),SIN))))));
```

> ##### 2. Prepare PAML input files.

```bash
#prepare seq file
mkdir gene
for i in ../align_NT_filter/*.fasta; do root=`basename $i .fasta`; echo "../scripts/paml_scripts/Fasta2Phylip.pl $i gene/${root}.phy";done > fa2phy.sh
ParaFly -c fa2phy.sh -CPU 30

#prepare parameter file
mkdir run_paml
perl ../scripts/paml_scripts/02.prepare.paml.Creat.CTL.pl gene run_paml

find `pwd` -name codeml.ctl > 1.txt; sed -e "s/^/cd &/g" 1.txt > paml.sh ; sed -i "s/codeml.ctl/; codeml ; cd ..\/..\/..\//g" paml.sh ; rm 1.txt;

#run paml
ParaFly -c paml.sh -CPU 30
```

> ##### 3. Extract paml results

```bash
perl ../scripts/paml_scripts/03.collect.paml.Freeratio.pl ./run_paml
perl ../scripts/paml_scripts/04.collect.paml.Oneratio_and_Tworatio.pl ./run_paml
perl ../scripts/paml_scripts/05.collect.paml.Branchsite.pl ./run_paml

#The results of positively selected genes are in the Branchsite.ps.txt file, and the results of rapidly evolving genes are in the Oneratio_and_Tworatio.txt file.
```
#### - Conserved Non-coding Elements Analysis

> ##### PART ONE: The pipeline for identifying CNEs.

> ##### Step 0: MafFilter

```bash
# Copy
mkdir 0.maf
for i in `seq 1 22` X; do cp /work1/fanhzh/project/gp_project/1.CNE/1.multiz/Chr$i/tmp.17.maf 0.maf/Chr$i.maf; done
# Filter the .maf file according to the score and reference
mkdir 1.maf_filter
for i in `seq 1 22` X; do echo "mafFilter -minScore=20000 -needComp=HSA -overlap  0.maf/Chr$i.maf > 1.maf_filter/Chr$i.filter.maf"; done > 01.filter.sh
perl easy_submit.pl -ss SLURM -t 1 -m 10G -l 1 -s yes -p 01.filter.sh
```

> ##### Step 1: Mafsort

```bash
# Sort the .maf file
mkdir 02.maf_sort
for file in 1.maf_filter/*.maf; do
    root=`basename $file .filter.maf`
    echo "maf-sort $file > 02.maf_sort/${root}.sort.maf"
done > 02.sort.sh
perl easy_submit.pl -ss SLURM -t 1 -m 10G -l 1 -s yes -p 02.sort.sh
```

> ##### Step 2: Get_ss

```bash
# Transfer the .maf format to .ss format
mkdir 03.maf2ss
for file in 02.maf_sort/*.maf; do
    root=`basename $file .sort.maf`
    echo "msa_view $file --in-format MAF --out-format SS > 03.maf2ss/${root}.ss"
done > 03.ss.sh
perl easy_submit.pl -ss SLURM -t 1 -m 10G -l 1 -s yes -p 03.ss.sh
```

> ##### Step 3: PhyloFit

```bash
# Rename the .gff file and split the file according to the chromosome number
sed -i "s/Chr/HSA.Chr/g" HSA.genome.gff
for i in `seq 1 22` X; do grep -w HSA.Chr$i HSA.genome.gff > HSA.Chr$i.gff; done

# Extract the 4DTV site for each chromosomes
mkdir 04.4dtv
for i in `seq 1 22` X
do
echo "msa_view 02.maf_sort/Chr${i}.sort.maf --4d --features HSA.Chr${i}.gff > 04.4dtv/4d-codons.Chr${i}.ss && msa_view 04.4dtv/4d-codons.Chr${i}.ss --in-format SS --out-format SS --tuple-size 1 > 04.4dtv/4d-sites.Chr${i}.ss"
done > 04.4dtv.sh
perl easy_submit.pl -ss SLURM -t 1 -m 10G -l 1 -s yes -p 04.4dtv.sh

# Combine all .ss file
msa_view --unordered-ss --out-format SS --aggregate HSA,BTA,MPE,FCA,PBE,AJU,LCA,PTI,NNE,PCR,CLU,VLA,UAR,UMA,UAM,AME,NVI,MER,LLU,SIN 04.4dtv/4d-sites.*.ss > 04.4dtv/all-4d.sites.ss

# Estimate the non-conserved model, the "nonconserved-all-4d" is equal to init.mod for subsequently analysis
phyloFit --tree "(HSA,(BTA,(MPE,((((((FCA,PBE),AJU),LCA),(PTI,NNE)),PCR),((CLU,VLA),((((UAR,UMA),UAM),AME),(((NVI,MER),LLU),SIN)))))))" --msa-format SS --out-root nonconserved-all-4d 04.4dtv/all-4d.sites.ss
```

> ##### Step 4: PhastCons (rate-limiting step)

```bash
# Estimate the conserved model and non-conserved model for each chromosomes
mkdir 05.PhastCons && cd 05.PhastCons
for file in ../03.maf2ss/*.ss ; do
    root=`basename $file .ss`
    echo "phastCons --target-coverage 0.3 --expected-length 45 --gc 0.4 --estimate-trees $root $file ../nonconserved-all-4d.mod --no-post-probs"
done > 05.phastCons.sh
perl easy_submit.pl -ss SLURM -t 1 -m 50G -l 1 -s yes -p 05.phastCons.sh
cd ..
```

> ##### Step 5: PhyloBoot

```bash
# Calculate the ave.cons.mod and ave.noncons.mod
mkdir 06.phyloboot && cd 06.phyloboot
ls ../05.PhastCons/*.cons.mod > cons.txt
phyloBoot --read-mods '*cons.txt' --output-average ave.cons.mod
ls ../05.PhastCons/*.noncons.mod > noncons.txt
phyloBoot --read-mods '*noncons.txt' --output-average ave.noncons.mod
```

> ##### Step 6: Get .bed file of conserved element 

```bash
mkdir elements scores
for file in ../03.maf2ss/*.ss; do
    root=`basename $file .ss`
    echo "phastCons --target-coverage 0.3 --expected-length 45 --most-conserved elements/${root}.bed --score $file ave.cons.mod,ave.noncons.mod > scores/${root}.wig"
done > get_bed.sh
perl easy_submit.pl -ss SLURM -t 1 -m 30G -l 1 -s yes -p get_bed.sh

# Modifying the format of the .bed file
cat elements/*.bed | sort -k1,1 -k2,2n > most-conserved.bed
awk '{OFS="\t"} $4="0",$5="0"' most-conserved.bed > replaced_most-conserved.bed
```

> ##### Step 7: Filter elements overlap the coding regions

```bash
gffread HSA.genome.gff -T -o HSA.genome.gtf

gtfToGenePred -genePredExt HSA.genome.gtf HSA.genome_refGene.txt

perl retrieve_seq_from_fasta.pl -format refGene -seqfile HSA.genome.fna ./HSA.genome_refGene.txt -out HSA.genome_refGeneMrna.fa

#The annotate_variation.pl script can be used to download humandb.tar.gz.
perl annotate_variation.pl -downdb -webfrom annovar refGene humandb/

cp HSA.genome_refGene.txt humandb/

perl annotate_variation.pl --outfile CNE -buildver HSA.genome replaced_most-conserved.bed humandb/

cp CNE.variant_function CNE.noExonandUTR.bed

sed -i "/exonic/d" CNE.noExonandUTR.bed
sed -i "/UTR3/d" CNE.noExonandUTR.bed
sed -i "/UTR5/d" CNE.noExonandUTR.bed
```

> ##### Step 8: Deal with CNEs bed

```bash
awk -F "\t" '{print $3"\t"$4"\t"$5}' CNE.noExonandUTR.bed > CNEs.ori.bed
# Merge beds closer than 10bp
bedtools merge -i CNEs.ori.bed -c 1 -d 10 -o count > CNEs.merge.bed
# Filter the CNEs length less than 30bp or 50bp
awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$3-$2}' CNEs.merge.bed > CNEs.cal.bed
awk '{if ($4>29) print $0}' CNEs.cal.bed > CNEs.filter.bed
awk '{if ($4>49) print $0}' CNEs.cal.bed > CNEs.50bp.bed
# Rename all the CNEs
awk 'BEGIN {FS="\t"; OFS="\t"} {$4="CNE"NR; print}' CNEs.filter.bed > CNEs.final.bed
# Statistics for CNEs length
cat CNEs.filter.bed | awk 'BEGIN{sum=0}{sum += $4}END{print sum}'
```

> ##### PART TWO: The pipeline for ForwardGenomics

```````bash
# ForwardGenomics
# create binary files for each genome
mkdir 0.2bit
ls /work1/fanhzh/project/gp_project/0.genome/simple/*.fa > genome.list
for i in `cat genome.list`; do root=`basename $i .genome.fa`;echo "mkdir 0.2bit/${root}; faToTwoBit $i 0.2bit/${root}/${root}.2bit";done > work.sh
perl easy_submit.pl -ss SLURM -t 1 -m 10G -l 1 -s yes -p work.sh

# deal maf file
mkdir 1.mafAddIRows
for i in /work1/fanhzh/project/gp_project/1.CNE/2.deal_maf/02.maf_sort/*.maf; do root=`basename $i .sort.maf`; echo "mafAddIRows $i 0.2bit/HSA/HSA.2bit 1.mafAddIRows/${root}.addI.maf";done > addI.sh

faSize -detailed HSA.genome.fna > HSA.chrom.sizes
for i in 1.mafAddIRows/*.maf; do root=`basename $i`;echo "~/Software/TOGA/CESAR2.0/tools/mafIndex $i 1.mafAddIRows/${root}.bb chromSizes=HSA.chrom.sizes";done > index.sh

#split maf file according to CNEs
mkdir 2.CNEs && cd 2.CNEs
for i in `seq 1 22` X; do grep -w "Chr$i" ../../06.phyloboot/CNEs.final.bed > Chr$i.CNEs.bed; done
for i in ./*.CNEs.bed; do root=`basename $i .CNEs.bed`; echo "/work1/fanhzh/project/gp_project/software/CESAR2.0-1.01/tools/mafExtract -regionList=$i -leaveEdgeMeta ../1.mafAddIRows/${root}.addI.maf.bb ${root}.CNEs.maf"; done > mafExtract.sh
perl easy_submit.pl -ss SLURM -t 2 -m 15G -s yes -p normal -l 1 mafExtract.sh

for i in ./*.CNEs.bed; do root=`basename $i .bed` ; awk '{print $4}' $i > ${root}.ID.txt; done

mkdir split_maf && cd split_maf
for i in ../*.CNEs.maf; do root=`basename $i .CNEs.maf`; echo "mkdir ${root}_split_maf; cd ${root}_split_maf; mafSplit ../../${root}.CNEs.bed ${root} ../../${root}.CNEs.maf; cd ../" ; done > mafSplit.sh
perl easy_submit.pl -ss SLURM -t 2 -m 15G -s yes -p normal -l 1 mafSplit.sh

for i in ./*split_maf; do rm $i/*.00.maf; done

# rename maf file
for i in ./*split_maf; do root=`basename $i _split_maf`; echo "bash generate_rename.sh $i ${root}.CNEs.ID.txt"; done > generate_rename_shell.sh

perl easy_submit.pl -ss SLURM -t 2 -m 15G -s yes -p normal -l 1 generate_rename_shell.sh

for i in ./*.batch_rename.sh; do echo "bash run_batch_rename.sh $i"; done > multi_run_batch_rename.sh

perl easy_submit.pl -ss SLURM -t 10 -m 30G -s yes -p normal -l 1 multi_run_batch_rename.sh

# Generate input files for ForwardGenomics
vi tree.nh
(HSA:0.0795224,((((PCR:0.0758239,(((AJU:0.00863126,LCA:0.00746509):0.000481099,(PBE:0.00747653,FCA:0.00705404):0.000944491):0.00319803,(NNE:0.00645053,PTI:0.00557254):0.00480762):0.0420401):0.0416664,((VLA:0.0132252,CLU:0.0114543):0.0826358,(((LLU:0.0243255,(MER:0.0163818,NVI:0.0168585):0.0142702):0.0502323,SIN:0.0774389):0.0199765,(AME:0.0203341,(UAM:0.00309861,(UAR:0.0018908,UMA:0.00227361):0.00155236):0.0162024):0.0421065):0.0202447):0.0196055):0.0472435,BTA:0.176567):0.00940529,MPE:0.133561):0.0795224);  #this tree in nonconserved-all-4d.mod
for i in `seq 1 22` X Y; do mkdir Chr${i}_prank; cd Chr${i}_prank; for j in ../Chr${i}_split_maf/*.maf; do root=`basename $j .maf`; echo "~/software/ForwardGenomics-1.0/scripts/Maf2SpanningSeq_PRANK.perl $j ${root} -runPrank -treeFile ../tree.nh -BDBFile Chr${i}.bdb -twoBitPath /work1/fanhzh/project/gp_project/1.CNE/2.deal_maf/07.ForwardGenomics/0.2bit"; done > Chr${i}.prank.sh; cd ../; done

#split scripts
for i in `seq 1 22` X; do cd Chr${i}_prank && split -l 20000 ./Chr${i}.prank.sh split_ -d && cd ..; done
for i in ./*prank; do mkdir $i/split; done
for i in `seq 1 22` X; do mv Chr${i}_prank/split_* Chr${i}_prank/split; done

for i in `seq 1 22` X; do for j in Chr${i}_prank/split/split_*; do root=`basename $j split_`; sed -i 's/\.\.\//\.\.\/\.\.\//g' $j; sed -i 's/-BDBFile Chr${i}.bdb/-BDBFile ${root}_Chr${i}.bdb/g' $j; done; done > replace_output.sh
bash replace_output.sh

# generate local and global files
tree_doctor -a tree.nh > tree.label
for i in `seq 1 22` X; do for j in Chr${i}_prank/split/split_[0-1][0-9]; do echo "sed -E 's|Maf2SpanningSeq_PRANK\.perl ../../Chr${i}_split_maf/(CNE[0-9]+)\.maf \1 .*-BDBFile (split_([0-9]+)_Chr${i}\.bdb).*|GetGlobalAndLocalPercentID.perl  \2 \1 -treeFile ../../tree.label  -allowedAncestralNodes HSA-PCR -requireNoOutgroup -local -global -GlobalBDBFile global_\3.bdb  -LocalBDBFile local_\3.bdb  -doNotIgnoreN|' $j > ${j}_gl"; done; done > replace.sh

bash replace.sh

for i in `seq 1 22` X; do for j in Chr${i}_prank/split/*gl; do echo "bash $j"; done; done > run_gl.sh
perl easy_submit.pl -ss SLURM -t 1 -m 30G -s yes -p normal -l 1 run_gl.sh

# deal Global and local files
ReadBDB.perl -brief global.bdb > global.txt
head -1 global.txt > global.header
sed -i "/species/d" global.txt
cat global.header global.txt > global.final.txt
ReadBDB.perl -brief local.bdb > local.txt
head -1 local.txt > local.header
sed -i "/branch/d" local.txt
cat local.header local.txt > local.final.txt

for i in `seq 1 22` X; do for j in Chr${i}_prank/split/global_[0-1][0-9].bdb; do  root=`basename $j .bdb`; echo "/public/home/huangchen/software/ForwardGenomics-1.0/scripts/ReadBDB.perl -brief $j > Chr${i}_prank/split/${root}.txt && head -1 Chr${i}_prank/split/${root}.txt > Chr${i}_prank/split/${root}.header && sed -i \"/species/d\" Chr${i}_prank/split/${root}.txt && cat Chr${i}_prank/split/${root}.header Chr${i}_prank/split/${root}.txt > Chr${i}_prank/split/${root}.final.txt"; done; done > deal_global.sh

bash deal_global.sh

for i in `seq 1 22` X; do for j in Chr${i}_prank/split/local_[0-1][0-9].bdb; do  root=`basename $j .bdb`; echo "/public/home/huangchen/software/ForwardGenomics-1.0/scripts/ReadBDB.perl -brief $j > Chr${i}_prank/split/${root}.txt && head -1 Chr${i}_prank/split/${root}.txt > Chr${i}_prank/split/${root}.header && sed -i \"/branch/d\" Chr${i}_prank/split/${root}.txt && cat Chr${i}_prank/split/${root}.header Chr${i}_prank/split/${root}.txt > Chr${i}_prank/split/${root}.final.txt"; done; done > deal_local.sh

bash deal_local.sh

#CNE element ID
for i in `seq 1 22` X; do for j in Chr${i}_prank/split/global_[0-1][0-9].final.txt; do root=`basename $j .final.txt`; tail -n +2 $j | cut -f1 -d " " > Chr${i}_prank/split/${root}.ID.txt; done; done

#The phenotype input file should be formatted with space-delimited binary values (1: phenotype present, 0: phenotype absent)
vi phenotype_AME.txt
species pheno
AJU 1
AME 0
BTA 1
CLU 1
FCA 1
HSA 1
LCA 1
LLU 1
MER 1
MPE 1
NNE 1
NVI 1
PBE 1
PCR 1
PTI 1
SIN 1
UAM 1
UAR 1
UMA 1
VLA 1


for i in `seq 1 22` X; do mkdir Chr${i}_prank/out; for j in Chr${i}_prank/split/global_[0-1][0-9].final.txt; do root=`basename $j .final.txt`; echo "/public/home/huangchen/software/ForwardGenomics-1.0/forwardGenomics.R --tree=./tree.label --elementIDs=Chr${i}_prank/split/${root}.ID.txt --listPheno=./phenotype_AME.txt --globalPid=${root}.final.txt --localPid=${root}.final.txt  --outFile=Chr${i}_prank/out/${root}.output.txt --method=all --weights=/public/home/huangchen/software/ForwardGenomics-1.0/lookUpData/branchWeights_CNE.txt --expectedPerIDs=/public/home/huangchen/software/ForwardGenomics-1.0/lookUpData/expPercentID_CNE.txt --minLosses=1 --thresholdConserved=0" | sed 's/--localPid=global/--localPid=local/'; done; done > run_ForwardGenomics.sh

perl easy_submit.pl -ss SLURM -t 1 -m 15G -p fat -s yes -l 1 run_ForwardGenomics.sh

cat Chr*_prank/out/*.txt > ForwardGenomics_out.txt

#Filter results with user-defined thresholds (adjust according to your needs)
awk 'NR == 1 {print; next};{gls = $7 + 0; pearson = $NF + 0};$3 == 19 && $5 == 37 && gls != NA && pearson != NA && pearson <=0.01' ForwardGenomics_out.txt > ForwardGenomics_out.filtered_0.01.txt

#ForwardGenomics shows high sensitivity, resulting in numerous CNE deletions in giant pandas. However, since some CNEs are also missing in other species, these require filtering.
cut -d ' ' -f 1 ForwardGenomics_out.filtered_0.01.txt > RECNEs.ID.txt
python filter_RECNEs.py  #need to revise parameters 
```````

> #### PART THREE: The pipeline for Phyloacc

> ##### Step 1: Prepare the input file

```bash
for i in `seq 1 22` X
do
mkdir Chr${i}
grep -w Chr${i} ./CNEs.final.bed > Chr${i}/Chr${i}.bed
cat << EOF >> Chr${i}/Chr${i}.sh
perl 01.convertMaf2List.pl /work1/fanhzh/project/gp_project/1.CNE/2.deal_maf/02.maf_sort/Chr${i}.sort.maf Chr$i    #need to revise parameters 
perl 02.lst2fas.pl Chr${i}.maf.lst ./Chr${i}.bed
EOF
done

for i in in `seq 1 22` X; do echo "cd $i && bash $i.sh && cd ..";done > run.sh
perl easy_submit.pl -ss SLURM -t 1 -m 15G -p fat -s yes -l 1 run.sh
```

> ##### Step 2: Split the large number of CNEs to several small directory

```bash
for i in `seq 1 22` X; do echo "find Chr$i/ -name \"*.fa\" | xargs -i cp {} CNEs"; done > move.sh     #Copy all chromosomal CNE FASTA files to the CNEs directory
perl easy_submit.pl -ss SLURM -t 1 -m 15G -p fat -s yes -l 1 move.sh
```

> ##### Step 3: Estimate the neutral mod for phyloacc

```bash
vi tree.txt
(HSA,(BTA,(MPE,((((((FCA,PBE),AJU),LCA),(PTI,NNE)),PCR),((CLU,VLA),((((UAR,UMA),UAM),AME),(((NVI,MER),LLU),SIN)))))));

phyloFit --tree tree.nwk --init-random --subst-mod SSREV --out-root neut_ver1_1 --msa-format SS --sym-freqs --log phyloFit_ver1_1.log ../04.4dtv/all-4d.sites.ss &> phyloFit_ver1_1.out

phyloFit --init-model neut_ver1_1.mod --out-root neut_ver1_1_update --msa-format SS --sym-freqs --log neut_ver1_1_update.log ../04.4dtv/all-4d.sites.ss &> neut_ver1_1_update.out

tree_doctor --name-ancestors neut_ver1_1_update.mod > neut_ver1_1_FINAL.mod
mv neut_ver1_1_FINAL.mod input.mod
```

> ##### Step 4: Run the phyloacc under the specific conda environment

```bash
conda activate phyloacc
mkdir phyloP_out
 for i in in `seq 1 22` X; do echo "phyloP --method LRT --features /work1/fanhzh/project/gp_project/1.CNE/2.deal_maf/07.ForwardGenomics/2.CNEs/Chr$i.CNEs.bed --mode CONACC --branch AME ./neut_ver1_1_FINAL.mod --msa-format MAF ../02.maf_sort/Chr$i.sort.maf > phyloP_out/Chr1.AME.out";done > phyloP.sh
 perl easy_submit.pl -ss SLURM -t 1 -m 15G -p fat -s yes -l 1 phyloP.sh
```
#### - Pseudogene Identification
> ##### 1. Preparation of input files
```
## At least four input files are required
## 1. query.2bit; 2. target.2bit; 3. all_clean.chain; 4.Species.gene.bed (bed12 format)
## 1, 2, 3 come from whole genome alignment pipeline

### Bed12
gff3ToGenePred HSA.genome.gff HSA.genePred
genePredToBed HSA.genePred HSA.bed
awk '{printf "sed -i \"s/"$1"/"$2"/g\" HSA.bed\n"}' HSA.seqname.list | sh
sed -i "/^NC\|^NW/d" HSA.bed
sort -k 1,1V -k 2,2n HSA.bed > HSA.sort.bed
for i in `seq 1 22` X Y; do grep -w "Chr$i" HSA.sort.bed > Chr${i}.sort.bed; done

### isoforms
awk '{if ($3 == "mRNA") print $0}' HSA.genome.gff > HSA.mRNA.bed
sed "s/;/ /g" HSA.mRNA.bed | awk '{print $10,$9}' - | sed "s/ID=\|Parent=gene-//g" - | perl -p -i -e 's/ /\t/g' > Isoforms.txt

### 2bit,chain
mkdir HSA_AME_chain HSA_CLU_chain HSA_FCA_chain HSA_MER_chain HSA_NSC_chain HSA_PCR_chain HSA_UMA_chain
for i in `seq 1 22` X Y; do 
cp /work1/fanhzh/project/gp_project/1.CNE/Chr$i/Chr$i_AME*/all_clean.chain HSA_AME/Chr$i.all.clean.chain;done   #Copy chain files of each species to the corresponding directories respectively
cp /work1/fanhzh/project/gp_project/1.CNE/2.deal_maf/07.ForwardGenomics/0.2bit/*/*.2bit ./
mkdir HSA_2bit
faSplit byName HSA.genome.fna HSA_2bit
for i in HSA_2bit/*.fa;do root=`basename $i fa`; faToTwoBit $i HSA_2bit/${root}2bit;done
```
> ##### Run TOGA

```bash
conda activate TOGA

# run TOGA
for i in `seq 1 22` X Y; do echo "python ~/software/TOGA/toga.py --isoforms Isoforms.txt --nextflow_config_dir ~/software/TOGA/nextflow_config_files --project_name AME_Chr$i ./HSA_AME_chain/Chr${i}_all_clean.chain ./Chr${i}.sort.bed.filtered.bed ./HSA_2bit/Chr${i}.2bit ./AME.2bit";done > AME_TOGA.sh
perl easy_submit.pl -ss SLURM -t 30 -m 200G -p long -l 1 -s yes AME_TOGA.sh

# filter results
awk '{if ($1 == "GENE" && $3 == "L" || $3 == "UL") print $0 }' AME_Chr*/loss_summ_data.tsv > AME_loss.txt

#extract specific loss in giant panda
grep -Fxv -f <(cat CLU_loss.txt  FCA_loss.txt  UMA_loss.txt  NSC_loss.txt MER_loss.txt | sort -u) AME_loss.txt > AME_uniq_loss.txt
```
#### - Functional Enrichment Analysis
  ```
  ```
