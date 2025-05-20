$sorg="Chr$ARGV[0]\.Sorghum_ancestral_allele_V3.1.txt";
open(T,"$sorg")||die "Can't open your file!";
while(<T>)
{
	chomp;
@aa=split;##5       367     T       1       G       0       T Chr     Pos     allele1 Freq1   allele2 Freq2   ancestral_allele
$mk="$aa[0]\_$aa[1]";
$hash{$mk}=$aa[6];
	}
close(T);

$infilename="Chr$ARGV[0]\.landrace_107samples.recode.vcf";
open(T,"$infilename")||die "Can't open your file!";
open(R,">Chr$ARGV[0]\.landrace_107samples.propinquum.vcf");
while(<T>)
{
	chomp;
	if(/^##/){print R "$_\n";next;}
	if(/CHROM/){
	@aa=split /\t/,$_,10;
	print R "$aa[0]	$aa[1]	$aa[2]	$aa[3]	$aa[4]	$aa[5]	$aa[6]	$aa[7]	$aa[8]	propinquum	$aa[9]\n";next;
	}
@aa=split /\t/,$_,10;
$mk1="$aa[0]\_$aa[1]";
if(exists($hash{$mk1}))
{
$g=$hash{$mk1};
if($aa[3] eq $g)
{
$ge="0|0";
print R "$aa[0]	$aa[1]	$aa[2]	$aa[3]	$aa[4]	$aa[5]	$aa[6]	$aa[7]	$aa[8]	$ge	$aa[9]\n";
}
if($aa[4] eq $g)
{
$ge="1|1";
print R "$aa[0]	$aa[1]	$aa[2]	$aa[3]	$aa[4]	$aa[5]	$aa[6]	$aa[7]	$aa[8]	$ge	$aa[9]\n";
}
}
}