use warnings;
use strict;

my $annot = shift or die "Usage:\nperl $0  RPM_annot  > mod_RPM_annot\n";
open(F,$annot) || die "Can't open the $annot!\n";
my $n = 0;
$annot =~ /-(\d+).fa/;
$n = $1;
#print STDERR "$n\n";
while(<F>){
    next if(/^pValue/);
    my @a = split /\s+/;
    if($n){
        $a[4] += $n;
        $a[5] += $n;
    }
    s/^(\S+\s+\S+\s+\S+\s+\S+\s+)\S+(\s+)\S+/$1$a[4]$2$a[5]/ if($n);
    print;
}
close F;
