#!/bin/bash

snakemake -p --profile slurm --conda-prefix $HOME/.conda/snakemake
