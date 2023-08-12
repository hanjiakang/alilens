#!/bin/bash
# Run CMB and noise simulations
python -W ignore sims.py -np 70

# Run reconstruction process
python plot.py

# Generate reconstruction potential map based on the klm
python klm.py

