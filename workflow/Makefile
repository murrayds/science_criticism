all:
	snakemake -j 4

download:
	snakemake -j 1 download --rerun-triggers mtime --use-conda --conda-frontend conda

novelty:
	snakemake -j 2 novelty --rerun-triggers mtime --use-conda --conda-frontend conda

novelty-hpc:
	snakemake -j 20 novelty --force --use-conda --conda-frontend conda

match:
	snakemake -j 2 match --rerun-triggers mtime --use-conda --conda-frontend conda

results:
	snakemake -j 4 results --rerun-triggers mtime --use-conda --conda-frontend conda
