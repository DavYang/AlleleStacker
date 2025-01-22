#!/bin/bash

# Source conda
source ~/.bashrc

# Activate conda environment
conda activate anc_vig

# Run the tests
PYTHONPATH=. python -m pytest test_variant_mapper_4.py -v
