#/bin/perl
use strict;
die "Usage:\nperl $0  <aligned_seq_dir(XXX.fasta)>  <outdir>" if(@ARGV<2);

my $aligned_seq_dir = shift(@ARGV);
my $outdir = shift(@ARGV);
my @order=@ARGV;
my $outfile = "merged.aln.fasta";

opendir(PRDIR, "$aligned_seq_dir") or die "$!"; my @files = grep(/fasta$/, readdir(PRDIR)); close PRDIR;
my %hash; my $id;
foreach my $file (sort (@files)){
	open(IN, "$aligned_seq_dir/$file") or die "$!";
	while(my $line = <IN>){
		chomp($line);
		#if($line =~ m/^>(\D\D\D)_/){
		if($line =~ m/^>([a-zA-Z]*)/){
			$id = $1;
		}else{
			if($line =~ m/^\s+/){print "$file\twrong\n";}
			$hash{$id} .= $line;
		}
	}
	close IN;
}

open(MER, ">$outdir/$outfile") or die "$!";
foreach my $sp (sort (keys %hash)){
	print "$sp\n";
	print MER ">$sp\n$hash{$sp}\n";
}
close MER;
close OUT;

my $trim_cut = "trim.$outfile";
#system("trimal -in $outdir/$outfile -out $outdir/$trim_cut -automated1");
system("trimal -in $outdir/$outfile -out $outdir/$trim_cut -nogaps");

