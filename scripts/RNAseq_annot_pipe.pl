#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

########## Created by Xin Du in 2022/02/12 ##########

my ($fq_dir,$transcript,$genome,$genome_guided,$trinity_path,$cdhit_path,$run_pasa,$pasa_path,$hisat2_path,$species,$max_intronlen,$threads,$merge_bam,$stringtie,$stringtie_path,$help,$fl_trans,$alignment_tools,$exec_type,$partition);

GetOptions(
    'fq_dir|fd=s' => \$fq_dir,
    'transcript|trans=s' => \$transcript,
    'alignment_tool|at=s' => \$alignment_tools,
    'genome|g=s' => \$genome,
    'exec_type|et=s' => \$exec_type,
    'genome_guided|gg' => \$genome_guided,
    'full_length|fl=s' => \$fl_trans,
    'trinity_path|tp=s' => \$trinity_path,
    'cdhit_path|cp=s' => \$cdhit_path,
    'run_pasa|r=s' => \$run_pasa,
    'pasa_path|pp=s' => \$pasa_path,
    'hisat2_path|hp=s' => \$hisat2_path,
    'species|s=s' => \$species,
    'max_intronlen|m=i' => \$max_intronlen,
    'threads|t=i' => \$threads,
    'stringtie_path|sp=s' => \$stringtie_path,
	'partition|p=s' => \$partition,
    'help|h' => \$help
);

my $help_info = "
Usage: perl $0 [options]
###       Necessary      ###
--fq_dir|-fd        str   The PATH of saving fastq files (file name: *_R1|2.fastq.gz).
--transcript|-trans file  The full-length or already assembled transcript data(fasta file). This is necessary when --fq_dir no value.
--full_length|-fl   file  The full length transcript fasta file. This is necessary when --fq_dir and --transcript no value.

###        Options       ###
--exec_type|-et      str   options: local, submit (default).
--alignment_tool|-at str   PASA alignment tools, options: blat,gmap,minimap2, default: blat,gmap
--full_length|-fl    file  Full length transcript fasta file, default: off
--run_pasa|-r        str   Whether run PASA(yes: default or no)
--gemome|-g          file  reference genome for running guide trinity or PASA
--species|-s         str   The species name of reference genome [target]
--max_intronlen|-m   int   max_intronlen [500000]
--genome_guided|-gg        whether execute genome_guided Trinity.
--threads|-t         int   The number of threads for running Trinity, hisat2 or other software [30]
--partition|-p       str   partition for submit task [NA]
--trinity_path|-tp   str   The Trinity full path [ENV]
--pasa_path|-pp      str   PASA home dir path [ENV]
--hisat2_path|-hp    str   hisat2 home dir path [ENV]
--cdhit_path         str   cd-hit home dir path [ENV]
--stringtie_path|-st str   stringtie executable path [ENV]
--help|-h                  Display help information.

###   All dependent software could be read from your \$ENV by default   ###
";
die $help_info if((!$fq_dir and !$transcript and !$fl_trans) or $help);

$alignment_tools = "blat,gmap" if(!$alignment_tools);
$species = "target" if(!$species);
$max_intronlen = 500000 if(!$max_intronlen);
$run_pasa = "yes" if(!$run_pasa);
$threads = 30 if(!$threads);
$exec_type = "submit" if(!$exec_type);
$partition = "NA" if(!$partition);

#$stringtie_path = "stringtie" if(!$stringtie_path);
my $submit_script;
if($exec_type eq "submit"){
    chomp($submit_script = `which easy_submit.pl`);
    die "easy_submit.pl can't be found in your ENV!\n" if(!$submit_script);
}
chomp($pasa_path = `which Launch_PASA_pipeline.pl 2>/dev/null`) if(!$pasa_path);
die "Launch_PASA_pipeline.pl can't be found in your ENV!\n" if(!$pasa_path);
chomp($pasa_path = `realpath $pasa_path`);
$pasa_path =~ s/Launch_PASA_pipeline.pl$//;
chomp($cdhit_path = `which cd-hit 2>/dev/null`) if(!$cdhit_path);
die "cd-hit isn't found in your ENV!\n" if(!$cdhit_path);
$cdhit_path =~ s/cd-hit$//i;
## check samtools
chomp(my $samtools = `which samtools 2>/dev/null`);
die "The samtools can't find in the your ENV!\n" if(!$samtools);

chomp(my $work_path = `pwd`);
$genome = "NA" if(!$genome);
print STDERR "The genome-guided RNAseq assembly will be run!\n" if(-f $genome and $genome_guided);
my %fq;
my $data_type = "raw";

if($fq_dir){
    # identifying and classifying fastq.gz in the data folder
    my @all_fq1 = <$fq_dir/*_R1.*.gz>;
    my @j = grep {/fq|fastq/} @all_fq1;
    die "No *_R1|2.fastq|fq.gz can't be found in the $fq_dir!\n" if(@j == 0);
    for my $fq1(@all_fq1){
        next unless($fq1 =~ /fq|fastq/);
        my $fq2 = $fq1;
        $fq2 =~ s/R1.(fastq|fq)/R2.$1/;
        $fq1 =~ /$fq_dir\/(\S+)_R1.(fastq|fq)/;
        my $tissue = $1;
        $fq{$tissue} = [$fq1,$fq2];
    }
}elsif($transcript){
    $data_type = "done";
    print STDERR "The full-length or already assembled transcript ($transcript) will be used!\n";
}else{
	$data_type = "done";
	print STDERR "The full-length transcript ($fl_trans) will be used!\n";
}

if(-f $genome and $genome_guided and $fq_dir){
    # check dependent software path
    chomp($hisat2_path = `which hisat2 2>/dev/null`) if(!$hisat2_path);
    die "The hisat2 can't find in the your ENV!\n" if(!$hisat2_path);
    chomp(my $samtools = `which samtools 2>/dev/null`);
    die "samtools isn't found in your ENV!\n" if(!$samtools);
    $hisat2_path =~ s/hisat2$//;
    chomp($trinity_path = `which Trinity 2>/dev/null`) if(!$trinity_path);
    die "Trinity isn't found in your ENV!\n" if(!$trinity_path);
    chomp($trinity_path = `realpath $trinity_path`);
    $trinity_path =~ s/Trinity$//;
    # build genome index
    mkdir "genome_index";
    `$hisat2_path/hisat2-build $genome $work_path/genome_index/$species 2>/dev/null` if(!-f ".build_index.done");
    `touch ".build_index.done"` if(!-f ".build_index.done");
    open(O1,">01.1.run_hiast2.sh");
	open(O2,">01.2.run_GG_trinity.sh");
    for my $tissue(sort keys %fq){
        my ($fq1,$fq2) = ($fq{$tissue}->[0],$fq{$tissue}->[1]);
        print O1 "$hisat2_path/hisat2 -x $work_path/genome_index/$species -p $threads --fr --min-intronlen 20 --max-intronlen $max_intronlen --dta -1 $fq1 -2 $fq2 2>$tissue\_hisat2.log | samtools view -@ $threads -bS - | samtools sort -@ $threads - -o $work_path/${species}_${tissue}_RNAseq.bam\n";
		print O2 "$trinity_path/Trinity --genome_guided_bam $work_path/${species}_${tissue}_RNAseq.bam --max_memory 60G --genome_guided_max_intron $max_intronlen --CPU $threads --output $work_path/trinity_${species}_${tissue}_GG --full_cleanup\n";
    }
    close O1;
	close O2;
    `bash 01.1.run_hiast2.sh` if($exec_type eq "local" and !-f ".01_1.done");
    `$submit_script 01.1.run_hiast2.sh -t $threads -m 80G -s yes -e no -l 1 -jp $species\_align -p $partition` if(!-f ".01_1.done" and $exec_type eq "submit");
	`touch ".01_1.done"`;
	`bash 01.2.run_GG_trinity.sh` if($exec_type eq "local" and !-f ".01_2.done");
	`$submit_script 01.2.run_GG_trinity.sh -t $threads -m 80G -s yes -e no -l 1 -jp $species\_align -p $partition` if(!-f ".01_2.done" and $exec_type eq "submit");
	`touch ".01_2.done"`;
    # merge bam and run stringtie
    my @all_bam = <$work_path/$species*.bam>;
    my $all_bam = join " ", @all_bam;
    `samtools merge -@ 24 -o merged_sorted.bam $all_bam` if(!-f ".merge_bam.done");
    `touch .merge_bam.done`;
    chomp($stringtie_path = `which stringtie 2>/dev/null`) if(!$stringtie_path);
    die "stringtie isn't found in your ENV!\n" if(!$stringtie_path);
    chomp(my $transdecoder_path = `which TransDecoder.LongOrfs`);
    die "TransDecoder.LongOrfs isn't found in your ENV!\n" if(!$transdecoder_path);
	chomp($transdecoder_path = `realpath $transdecoder_path`);
    $transdecoder_path =~ s/TransDecoder.LongOrfs$//;
    `$stringtie_path merged_sorted.bam -o merged_sorted.bam.gtf -p $threads` if(!-f ".bam2gtf.done");
    `touch .bam2gtf.done`;
	if(!-f ".trans2cds.done"){
	    `$transdecoder_path/util/gtf_genome_to_cdna_fasta.pl merged_sorted.bam.gtf $genome > stringtie_transcripts.fasta`;
    	`$transdecoder_path/util/gtf_to_alignment_gff3.pl merged_sorted.bam.gtf > stringtie_RNAseq_pasa.gff3`;
	    `$pasa_path/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta stringtie_transcripts.fasta --pasa_transcripts_gff3 stringtie_RNAseq_pasa.gff3`;
	}
    `touch .trans2cds.done`;
    # filter and merge transcripts.fasta
    open(O,">02.GG_filter.sh");
    for my $tissue(sort keys %fq){
        my ($fq1,$fq2) = ($fq{$tissue}->[0],$fq{$tissue}->[1]);
        print O "$trinity_path/util/misc/get_longest_isoform_seq_per_trinity_gene.pl $work_path/trinity_$species\_$tissue\_GG.Trinity-GG.fasta > $work_path/trinity_$species\_$tissue\_GG.unigene1.fasta; $cdhit_path/cd-hit-est -i $work_path/trinity_$species\_$tissue\_GG.unigene1.fasta -o $work_path/trinity_$species\_$tissue\_GG.cdhit.fasta -T $threads -M 10000 -c 0.9 -n 9 -d 0 ; perl ~dux21/tools/RNAseq_anno/rename_seq.pl $work_path/trinity_$species\_$tissue\_GG.cdhit.fasta GG_$tissue 6 > $work_path/trinity_$species\_$tissue\_GG.rename.fasta\n";
    }
    close O;
    `$submit_script 02.GG_filter.sh -t 5 -m 20G -s yes -e no -l 1 -jp filter -p $partition` if(!-f ".02_GG_filter.done" and $exec_type eq "submit");
    `bash 02.GG_filter.sh` if(!-f ".02_GG_filter.done" and $exec_type eq "local");
    `touch ".02_GG_filter.done"`;
    `cat $work_path/trinity_*_GG.rename.fasta > raw_GG_transcripts.fasta`;
	my $cd_hit_cmd = "$cdhit_path/cd-hit-est -i raw_GG_transcripts.fasta -o GG_transcripts.fasta -c 0.9 -n 9 -d 0 -M 30000 -T $threads";
	open(T,">02.GG_cdhit.sh");
	print T "$cd_hit_cmd\n";
	close T;
    `$submit_script 02.GG_cdhit.sh -t $threads -m 40G -s yes -e no -jp cdhit` if(!-f ".GG_cdhit.done");
    `touch ".GG_cdhit.done"`;
    &denovo_asm();
    `cat GG_transcripts.fasta denovo_final_transcripts.fasta > All_transcripts.fasta` if(!$fl_trans);
    if($fl_trans){
        `cat $fl_trans GG_transcripts.fasta denovo_final_transcripts.fasta > All_transcripts.fasta`;
        `perl -alne 'print \$1 if(/^>(\\S+)/);' $fl_trans > FL_cDNA.name`;
    }
    # run PASA
    exit 0 if($run_pasa eq "no");
    `mkdir -p "database"`;
    my $database = $work_path."/database/database.sqlite";
    #`cp $pasa_path/pasa_conf/pasa.alignAssembly.Template.txt pasa_align.config`;
    `perl -le 'for(\`cat \$ARGV[0]\`){chomp;s/<__DATABASE__>/\$ARGV[1]/;print}' $pasa_path/pasa_conf/pasa.alignAssembly.Template.txt $database > pasa_align.config`;
    open(O, ">03.run_GG_pasa.sh");
    print O "$pasa_path/Launch_PASA_pipeline.pl -c pasa_align.config -C -R -g $genome -t All_transcripts.fasta --TDN denovo_trans_name.list --trans_gtf merged_sorted.bam.gtf --ALIGNERS blat,gmap --CPU $threads -I $max_intronlen\n" if(!$fl_trans);
    print O "$pasa_path/Launch_PASA_pipeline.pl -c pasa_align.config -C -R -g $genome -t All_transcripts.fasta --TDN denovo_trans_name.list -f FL_cDNA.name --trans_gtf merged_sorted.bam.gtf --ALIGNERS $alignment_tools --CPU $threads -I $max_intronlen\n" if($fl_trans);
    close O;
    `$submit_script 03.run_GG_pasa.sh -t $threads -m 150G -s yes -e no -jp $species\_pasa -p $partition` if(!-f ".03_GG_PASA.done" and $exec_type eq "submit");
    `bash 03.run_GG_pasa.sh` if(!-f ".03_GG_PASA.done" and $exec_type eq "local");
    `touch ".03_GG_PASA.done"`;
	print STDERR "All completed, please check results!\n";
}

sub denovo_asm{
    # check dependent software path
    chomp($trinity_path = `which Trinity 2>/dev/null`) if(!$trinity_path);
    die "Trinity isn't found in your ENV!\n" if(!$trinity_path);
    chomp($trinity_path = `realpath $trinity_path`);
    $trinity_path =~ s/Trinity$//;
    chomp($cdhit_path = `which cd-hit 2>/dev/null`) if(!$cdhit_path);
    die "cd-hit isn't found in your ENV!\n" if(!$cdhit_path);
    $cdhit_path =~ s/cd-hit$//i;
    chomp($pasa_path = `which Launch_PASA_pipeline.pl 2>/dev/null`) if(!$pasa_path);
    die "Launch_PASA_pipeline.pl can't be found in your ENV!\n" if(!$pasa_path);
    $pasa_path =~ s/Launch_PASA_pipeline.pl$//;
    open(O,">01.run_denovo_trinity.sh");
    for my $tissue (sort keys %fq){
        my ($fq1,$fq2) = ($fq{$tissue}->[0],$fq{$tissue}->[1]);
        print O "$trinity_path/Trinity --seqType fq --max_memory 60G --left $fq1 --right $fq2 --CPU $threads --output $work_path/trinity_$species\_$tissue\_denovo --full_cleanup\n";
    }
    close O;
    my $t = $threads;
    `$submit_script 01.run_denovo_trinity.sh -t $t -m 60G -s yes -l 1 -e no -p $partition` if($exec_type eq "submit" and !-f ".01_denovo.done");
    `bash 01.run_denovo_trinity.sh` if($exec_type eq "local" and !-f ".01_denovo.done");
    `touch ".01_denovo.done"`;
    # merge transcripts.fasta
    open(O,">02.denovo_filter.sh");
    for my $tissue(sort keys %fq){
        my ($fq1,$fq2) = ($fq{$tissue}->[0],$fq{$tissue}->[1]);
        print O "$trinity_path/util/misc/get_longest_isoform_seq_per_trinity_gene.pl $work_path/trinity_$species\_$tissue\_denovo.Trinity.fasta > $work_path/trinity_$species\_$tissue\_denovo.unigene1.fasta; $cdhit_path/cd-hit-est -i $work_path/trinity_$species\_$tissue\_denovo.unigene1.fasta -o $work_path/trinity_$species\_$tissue\_denovo.cdhit.fasta -T $threads -M 30000 -c 0.9 -n 9 -d 0 ; perl ~dux21/tools/RNAseq_anno/rename_seq.pl $work_path/trinity_$species\_$tissue\_denovo.cdhit.fasta $tissue 6 > $work_path/trinity_$species\_$tissue\_denovo.rename.fasta\n";
    }
    close O;
    `$submit_script 02.denovo_filter.sh -t $threads -m 50G -s yes -e no -l 1 -jp filter -p $partition` if(!-f ".02_denovo.done" and $exec_type eq "submit");
    `bash 02.denovo_filter.sh` if(!-f ".02_denovo.done" and $exec_type eq "local");
    `touch ".02_denovo.done"`;
    if(-f "denovo_raw_transcripts.fasta"){
        chomp(my $seq_num = `grep -c ">" "denovo_raw_transcripts.fasta"`);
        `cat trinity_*_denovo.rename.fasta > denovo_raw_transcripts.fasta` if($seq_num == 0);
        chomp($seq_num = `grep -c ">" denovo_raw_transcripts.fasta`);
        die "No sequence in the denovo_raw_transcripts.fasta!\n" if($seq_num == 0);
    }else{
        `cat trinity_*_denovo.rename.fasta > denovo_raw_transcripts.fasta`;
    }
	open(T,">02.denovo_cdhit.sh");
	print T "$cdhit_path/cd-hit-est -i denovo_raw_transcripts.fasta -o denovo_final_transcripts.fasta -c 0.9 -n 9 -d 0 -M 30000 -T $threads\n";
	close T;
    `$submit_script 02.denovo_cdhit.sh -t $threads -m 40G -s yes -e no -jp cdhit` if(!-f ".denovo_cdhit.done");
    `perl -alne 'print "\$1" if(/^>(\\S+)/)' denovo_final_transcripts.fasta > denovo_trans_name.list`;
    `touch .denovo_cdhit.done`;
}

if(!$genome_guided and $fq_dir){
    &denovo_asm();
    if($run_pasa eq "yes"){
        next if($run_pasa eq "no");
        `mkdir -p "database"`;
        my $database = $work_path."/database/database.sqlite";
        #$database =~ s/\//\\\//g;
        `perl -le 'for(\`cat \$ARGV[0]\`){chomp;s/<__DATABASE__>/\$ARGV[1]/;print}' $pasa_path/pasa_conf/pasa.alignAssembly.Template.txt $database > pasa_align.config`;
        open(O,">03.run_denovo_pasa.sh");
        print O "$pasa_path/Launch_PASA_pipeline.pl -c pasa_align.config -C -R -g $genome -t denovo_final_transcripts.fasta --TDN denovo_trans_name.list --ALIGNERS $alignment_tools --CPU $threads -I $max_intronlen --TRANSDECODER\n";
        close O;
        `$submit_script 03.run_denovo_pasa.sh -t $threads -m 150G -s yes -e no -jp $species\_pasa -p $partition` if(!-f ".03_denovo.done" and $exec_type eq "submit");
        `bash 03.run_denovo_pasa.sh` if($exec_type eq "local" and !-f ".03_denovo.done");
        `touch ".03_denovo.done"`;
    }else{
        print STDERR "The Denovo transcript assembly have done!\n";
    }
}

if($transcript and !$fq_dir){
	die "$transcript isn't exists!\n" if(!-f $transcript);
	`mkdir -p "database"`;
	my $database = $work_path."/database/database.sqlite";
	`perl -le 'for(\`cat \$ARGV[0]\`){chomp;s/<__DATABASE__>/\$ARGV[1]/;print}' $pasa_path/pasa_conf/pasa.alignAssembly.Template.txt $database > pasa_align.config`;
	`perl -alne 'print \$1 if(/^>(\\S+)/)' $transcript > custom_trans.name`;
	if(!$fl_trans){
	    open(O, ">run_pasa.sh");
    	print O "$pasa_path/Launch_PASA_pipeline.pl -c pasa_align.config -C -R -g $genome -t $transcript --ALIGNERS $alignment_tools --CPU $threads -I $max_intronlen --TDN custom_trans.name --TRANSDECODER\n";
    	close O;
	}else{
		`perl -alne 'print \$1 if(/^>(\\S+)/)' $fl_trans > FL_cDNA.name`;
		`cat $transcript $fl_trans > All_transcripts.fasta`;
		open(O, ">run_pasa.sh");
		print O "$pasa_path/Launch_PASA_pipeline.pl -c pasa_align.config -C -R -g $genome -t All_transcripts.fasta --ALIGNERS $alignment_tools --CPU $threads -I $max_intronlen -f FL_cDNA.name --TDN custom_trans.name --TRANSDECODER\n";
		close O;
	}
	`$submit_script run_pasa.sh -t $threads -m 150G -s yes -e no -jp $species\_pasa -p $partition` if(!-f ".pasa.done" and $exec_type eq "submit");
	`bash run_pasa.sh` if(!-f ".pasa.done" and $exec_type eq "local");
	`touch ".pasa.done"`;
}

if($fl_trans and !$transcript and !$fq_dir){
	print STDOUT "You only give full-length RNAseq file, the PASA will be executed!\n";
	`mkdir -p "database"`;
	my $database = $work_path."/database/database.sqlite";
	`perl -le 'for(\`cat \$ARGV[0]\`){chomp;s/<__DATABASE__>/\$ARGV[1]/;print}' $pasa_path/pasa_conf/pasa.alignAssembly.Template.txt $database > pasa_align.config`;
	`perl -alne 'print \$1 if(/^>(\\S+)/)' $fl_trans > FL_cDNA.name`;
	open(O, ">run_pasa.sh");
	print O "$pasa_path/Launch_PASA_pipeline.pl -c pasa_align.config -C -R -g $genome -t $fl_trans -f FL_cDNA.name --ALIGNERS $alignment_tools --CPU $threads -I $max_intronlen\n";
	close O;
	`$submit_script run_pasa.sh -t $threads -m 150G -s yes -e no -jp $species\_pasa -p $partition` if(!-f ".pasa.done" and $exec_type eq "submit");
	`bash run_pasa.sh` if(!-f ".pasa.done" and $exec_type eq "local");
	`touch ".pasa.done"`;
}
