#!/usr/bin/perl

use strict;
use warnings;

# 读取bed文件
my $bed_file = "$ARGV[0]";
open(my $bed_fh, "<", $bed_file) or die "Cannot open bed file: $bed_file\n";

my @bed_data;
while (my $line = <$bed_fh>) {
    chomp $line;
    my @fields = split("\t", $line);
    push @bed_data, \@fields;
}
close $bed_fh;

# 读取fasta文件
my $fasta_file = "$ARGV[1]";
open(my $fasta_fh, "<", $fasta_file) or die "Cannot open fasta file: $fasta_file\n";

my %fasta_data;
my $current_header;
while (my $line = <$fasta_fh>) {
    chomp $line;
    if ($line =~ /^>(.+)/) {
        $current_header = $1;
        $fasta_data{$current_header} = "";
    } else {
        $fasta_data{$current_header} .= $line;
    }
}
close $fasta_fh;

# 提取序列并合并
my %merged_sequences;
foreach my $entry (@bed_data) {
    my ($chrom, $start, $end, $flag, $strand) = @{$entry};
    my $sequence = substr($fasta_data{$chrom}, $start-1, $end-$start+1);

    if ($strand eq '-') {
        $sequence = reverse_complement($sequence);
    }

    if (exists $merged_sequences{$flag}) {
        $merged_sequences{$flag} .= $sequence;
    } else {
        $merged_sequences{$flag} = $sequence;
    }
}

# 写入新的fasta文件
my $output_file = "$ARGV[2]";
open(my $output_fh, ">", $output_file) or die "Cannot open output file: $output_file\n";

foreach my $flag (keys %merged_sequences) {
    my $sequence = $merged_sequences{$flag};
    print $output_fh ">$flag\n$sequence\n";
}
close $output_fh;

sub reverse_complement {
    my ($sequence) = @_;
    $sequence = reverse $sequence;
    $sequence =~ tr/ACGTacgt/TGCAtgca/;
    return $sequence;
}

# Usage: perl Construct_Ancestor_Genome.pl Input.bed Species.fasta Ancestor.fasta