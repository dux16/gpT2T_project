use warnings;
use strict;

my $in = shift or die "Usage:\nperl $0 trf_dat_out > trf_dat_out.mod\n";
$in =~ /-(\d+).fa/;
my $n;
$n = $1 if($1);
if($n){
    if($n > 0){
        for(`cat $in`){
            next if($_ !~ /^(Sequence|Parameters|^\d+\s+\d+)/i);
            print if($_ !~ /^(\d+)\s+(\d+)/);
            if(/^(\d+)\s+(\d+)/){
                my ($s,$e) = ($1,$2);
                $s += $n;
                $e += $n;
                s/^\d+\s+\d+/$s $e/;
                print $_;
            }
        }
    }else{
        for(`cat $in`){
            next if($_ !~ /^(Sequence|Parameters|^\d+\s+\d+)/i);
            print;
        }
    }
}else{
    for(`cat $in`){
        next if($_ !~ /^(Sequence|Parameters|^\d+\s+\d+)/i);
        print;
    }
}
