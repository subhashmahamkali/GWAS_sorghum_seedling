open(E,">run.txt");
for($i=1;$i<=10;$i++)
{
 $out="filter_Chr$i\.slurm";
 open(R,">$out"); #
$cmd1="vcftools --vcf landrace_107samples.recode.vcf --chr $i --recode --out Chr$i\.landrace_107samples";
print R "#!/bin/sh
#SBATCH --job-name=filter_Chr$i
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=120:00:00
#SBATCH --mem=10G
#SBATCH --error=filter_Chr$i\.err
#SBATCH --output=filter_Chr$i\.out
#SBATCH --licenses=common
$cmd1
";
print E "sbatch $out\n";  

}