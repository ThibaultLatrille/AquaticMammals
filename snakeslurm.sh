#!/usr/bin/env bash
snakemake --rerun-incomplete --forcerun -k --printshellcmds -j 512 --cluster "sbatch -J Hyphy{params.name} -p cpu -N 1 -o ./slurm/%x.%j.out -e ./slurm/%x.%j.err --cpus-per-task=1 --mem={params.mem} -t {params.time}"
snakemake -s workflow/Relax.smk --rerun-incomplete --forcerun -k --printshellcmds -j 512 --cluster "sbatch -J Hyphy{params.name} -p cpu -N 1 -o ./slurm/%x.%j.out -e ./slurm/%x.%j.err --cpus-per-task=1 --mem={params.mem} -t {params.time}"
