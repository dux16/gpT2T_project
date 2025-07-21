#!/usr/bin/env perl
use warnings;
use strict;

my $phase_out = shift or die "Usage:\nperl $0 phase_out\n";

open(F,$phase_out) || die "Can't open the $phase_out!\n";

open(Mat,">maternal_new.cut");
open(Pat,">paternal_new.cut");
open(Amb,">ambiguous_new.cut");
while(<F>){
    chomp;
    next if(/^#|^\s*$/);
    my @a = split /\s+/;
    my $reads = $a[1];
    if($a[-2] > 3*$a[-1] or ($a[-2]>$a[-1] and $a[-2] > 30)){
        print Pat "$reads\n";
    }elsif($a[-1] > 3*$a[-2] or ($a[-1]>$a[-2] and $a[-1] > 30)){
        print Mat "$reads\n";
    }else{
        print Amb "$reads\n";
    }
}
close F;
close Mat;
close Pat;
close Amb;
