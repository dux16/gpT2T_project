#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);
use Getopt::Long;
 
my $help_requested;
my $file;
my $out="tau.txt";
my $r=3;
 
sub usage {
        print "Usage: $0 [options]\n";
        print "Options:\n";
        print "  -h, --help        Display this help message\n";
        print "  -i, --in <file> Specify an input file\n";
        print "  -o, --out <file> Specify an output file      default: tau.txt\n";
        print "  -r, --replicates <int> number of replicates    default: 3 \n";
}
if (@ARGV == 0) {
        usage();
        exit;
}
 
GetOptions(
        'h|help' => \$help_requested,
        'i|in=s' => \$file,
        'o|out=s' => \$out,
        'r|replicates=i' => \$r,
);
 
if ($help_requested or !$file) {
        usage();
        exit;
}

die "Can't open the cal_maleTau_TPM_tissue_name.txt!\n" if(!-f "cal_maleTau_TPM_tissue_name.txt");
chomp(my @ts_name = `cat cal_maleTau_TPM_tissue_name.txt`);
open FL, "$file" or die "cannot open the file $file\n";
open OUT, ">$out" or die "no output file name\n";
while(my $line = <FL>){
        chomp $line;
        my @array = split /\s+/,$line;
        my $gene = shift @array;
        my $n = scalar @array;
        my $tau;
        unless ($n % $r == 0 ){
            die "not be  divisible\n";
        }
        unless ($n > $r){
            die "Total expression can not more than or equal the number of eplicates ";
        }
        my @tpm_values;
        my @groups;
        foreach my $element (@array){
                push @groups, $element;
                if (@groups == $r) {
                        my $average = sum(@groups) / $r; 
                        push @tpm_values, $average;
                        @groups = ();
                }
        }
		# find max expression tissue
		my $max_index = 0;
		for my $i(0..$#tpm_values){
			$max_index = $i if($tpm_values[$i] > $tpm_values[$max_index]);
		}
        my $max_expression = max(@tpm_values);
		next if($max_expression <= 0);
        for my $tpm (@tpm_values){
				#print "$tpm\n";
                my $xi = $tpm / $max_expression;
                $tau += (1 - $xi);
        }
        $tau = $tau / (@tpm_values - 1);
        print OUT "$gene\t$tau\t$ts_name[$max_index]:$max_expression\n";
}
close FL;
close OUT;
sub max {
    return (sort {$a <=> $b} @_)[-1];
}
