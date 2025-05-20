open(T,"$ARGV[0]")||die "Can't open your file!";##Chr1.landrace_107samples.propinquum.vcf
$out=$ARGV[0];
$out=~s/vcf/input\.txt/;
open(R,">$out");
print R "physPos\tgenPos\tx\tn\n";
while(<T>)
{
	chomp;
	next if(/^#/);
	if(/\.\|\./){next;}
@aa=split /\t/,$_,11;
@bb=split /\t/,$aa[10];
$anc=substr($aa[9],0,1);
$x=0;
$n=0;
foreach $g(@bb)
{

($a1, $a2) = ($1, $2) if ($g =~ /^(\S)\|(\S)/);
$n += 2;
        if ($a1 ne $anc) {$x++;}
        if ($a2 ne $anc) {$x++;}
}
if($x==0){next;}
$gp=$aa[1]/1000000;
print R "$aa[1]	$gp	$x	$n\n";
	}