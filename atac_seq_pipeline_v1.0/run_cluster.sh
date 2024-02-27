
module load anaconda
source activate snakemake_env
snakemake -np -s combined.py --unlock
snakemake -c8 -s combined.py --cluster "sbatch -p cpuPartition -N 1 --mem=10000 --ntasks=1  --cpus-per-task=1 --job-name={rule} --output={rule}.%j.out --error={rule}.%j.err" --jobs 8 --rerun-incomplete

