open(E,">run.txt");
for($i=1;$i<=10;$i++)
{
	$out="put_propinqaum_in_vcf_Chr$i\.slurm";
	open(R,">$out"); #
$cmd="perl merge_propinqaum_genotype.pl $i";
print R "#!/bin/sh
#SBATCH --job-name=merge_propinqaum_genotype_Chr$i
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=120:00:00
#SBATCH --mem=20G
#SBATCH --error=merge_propinqaum_genotype_Chr$i\.err
#SBATCH --output=merge_propinqaum_genotype_Chr$i\.out
#SBATCH --licenses=common
$cmd";
print E "sbatch $out\n";		

}
