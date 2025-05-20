open(E,">run1.txt");
open(T,"vcf_files.txt")||die "Can't open your file!";
while(<T>)
{
	chomp;
	$o="$_\.slurm";
	open(R,">$o");
$cmd="perl get_input.pl $_";
print R "#!/bin/sh
#SBATCH --job-name=$_
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=120:00:00
#SBATCH --mem=10G
#SBATCH --error=$_\.err
#SBATCH --output=$_\.out
#SBATCH --licenses=common
$cmd";
print E "sbatch $o\n";
	}
close(T,R);	