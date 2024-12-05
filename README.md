# Code for the paper ``The origin, consequence, and visibility of criticism in science''

This repository contains code to reproduce all figures and tables used for the paper titled [``The origin, consequence, and visibility of criticism in science''](https://arxiv.org/abs/2412.02809).

To run the code, you will need to complete the following steps:

1. Copy the file `workflow/config.template.yaml` to a new file named `config.yaml`. In the new file, fill in the missing parameters, particularly `data_dir` (where to find and save data) and `fig_dir` (where to save figures). These should be paths on your local system.
    * Note that `data_dir` should already contain the folders `raw` and `downloaded` which will be accessed by the workflow. If `downloaded` does not exist, the workflow will attempt to query BigQuery for the data.
    * If you are are running this repository with access to CCNR-Success BigQuery data, then you will need to fill in the appropriate paths to the relevant data tables in `config.yaml`. Otherwise leave these blank.
    * Additionally, if you have access to CCNR-Success BigQuery, then you will need to ensure that the Google Cloud CLI has been installed and initialized on your local system (see [this link](https://cloud.google.com/sdk/docs/install) for instructions).
3. If you do not already have snakemake (>8.0.0) installed, then from the CLI in root directory, call `make create_conda_env` which create a conda environment with snakemake. Activate this environment before proceeding.
4. In the CLI, navigate to the `workflow/` directory. To execute the code, run `make results`, which will install the necessary conda environments and execute the snakemake workflow, processing any data necessary to create the result figures. The output will be saved in the `fig_path` that was set in the `config.yml` file.
    * Note that is you are running on an apple computer with an M1 chip, you may encounter an issue when the workflow attempts to install `r-matchit`, which is not currently built for the arm64 architecture. The solutions to this include look at each .yaml file in the `workflow/envs` directory and install these into a single conda environment, and through R use the `install.packages("matchit")` command to install it yourself through the CRAN package management system. Once the new environment is activated, use the command `snakemake -j 4 results --rerun-triggers mtime` in the `workflow/` directory to execure the code.
    * If the workflow attempts to download files from BigQuery (which will raise an error unless you are setup with access), or if it attepts to execute the novelty subworkflow (which will take a very very long time), then this means that the files are out of date and snakemake feels the need to update them. To resolve this, use the command `touch <path-to-project>/science_criticism/data/downloaded/*` and `touch <path-to-project>/science_criticism/data/derived/matched/*`, which will update their timestamps. Then re-run the workflow.


## Citation

```bibtex
@misc{chen2024criticism,
  Author = {Bingsheng Chen and Dakota Murray and Yixuan Liu and Albert-László Barabási},
  Title = {The origin, consequence, and visibility of criticism in science},
  Year = {2024},
  Eprint = {arXiv:2412.02809},
}
```
