all:
	snakemake -j 4

download:
	snakemake -j 1 download --rerun-triggers mtime

novelty:
	snakemake -j 2 novelty --rerun-triggers mtime

novelty-hpc:
	snakemake -j 20 novelty --force

match:
	snakemake -j 2 match --rerun-triggers mtime

results:
	snakemake -j 4 results --rerun-triggers mtime
