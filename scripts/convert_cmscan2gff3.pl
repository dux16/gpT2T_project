#!/usr/bin/perl

use strict;
use warnings;
# use Smart::Comments;

my $dirfile = shift; #文件夹路径
my $prefix = shift; #物种名字
die "Usage:\nperl $0 dirfile prefix\n" if(!$prefix or !$dirfile);

my %info = &read_seed();
print "Read seed file finished!\n";

open OUT, ">$prefix.rfam.gff";
my $mark = "000001";
my @tblout = <$dirfile/*.tblout>;
for my $file (@tblout){
    $file =~ /$dirfile\/(\S+).tblout/;
    my $chr = $1;
    #print $chr,"\n";
    open IN, $file || die "$!\n";
    while(<IN>){
    	chomp;
    	next if /^#/;
	    my @a = split/\s+/;
	    #next if ($info{$a[2]}{type} !~ /Gene;\s+(.*?);/);
        my $type;
        if(exists $info{$a[2]}{type} and $info{$a[2]}{type} =~ /Gene;\s+(.*?);/){
            $type = $1;
        }else{
            next;
        }
	    my $start = 0;
        my $end = 0;
    	if ($a[11] =~ /\+/){
	        $start = $a[9];
	        $end = $a[10];
    	}else {
	        $start = $a[10];
	        $end = $a[9];
	    }
    	my $gene_id = $prefix."_".$type."_".$mark;
	    print OUT "$chr\tcmsearch\t$type\t$start\t$end\t$a[-11]\t$a[11]\t.\tID=$gene_id;type=$info{$a[2]}{type};annotation=$info{$a[2]}{description};\n";
    	$mark ++;
    }
    close IN;
}
close OUT;


sub read_seed {
    my $reffile = "/public/home/duxin/software/02.genome_annotation/01.repetitive_elements_annotation/infernal-1.1.4-linux-intel-gcc/Rfam_lib/Rfam.seed";
    my %rfam;
    open(REF, $reffile) || die "$!\n";
    my $rfam_id;
    while(<REF>){
    	chomp;
	    if(/\#=GF\s+AC\s+(\S+)/){
	        $rfam_id = $1;
    	    # print "$rfam_id";exit;
	    }
    	if (/\#=GF\s+ID\s+(.+)/){
	        $rfam{$rfam_id}{name} = $1;
	    }
    	if(/\#=GF\s+TP\s+(.+)/){
	        $rfam{$rfam_id}{type} = $1;
	    }
    	$rfam{$rfam_id}{description} = $1 if(/\#=GF\s+DE\s+(.+)/);
    }
    close REF;
    return %rfam;
}
