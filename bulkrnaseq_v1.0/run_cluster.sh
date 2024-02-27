
source activate snakemake_env
snakemake -c4 -s combined.py --cluster "sbatch -p amd_512 -N 1 --mem=4000 --ntasks=1  --cpus-per-task=2 --job-name={rule} --output={rule}.%j.out --error={rule}.%j.err" --jobs 2

