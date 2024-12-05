PROJ_NAME=criticism

.PHONY: reate_conda_env build_conda_env report create_env create_ipykernel

create_conda_env:
	conda create -n $(PROJ_NAME) -c bioconda -c conda-forge python=3.11 snakemake=8.16.0

export_conda_env:
	conda env export > environment.yml

build_conda_env:
	conda env create -f environment.yml

create_env:
	python3 -m venv $(PROJ_NAME)-env

create_ipykernel:
	python3 -m ipykernel install --user --name=$(PROJ_NAME)
