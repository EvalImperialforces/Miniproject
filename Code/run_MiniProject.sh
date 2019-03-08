#!/bin/bash
# Author: Eva Linehan el1718@ic.ac.uk
# Script: Run Miniproject
# Desc: Script to run Miniproject from code to final LaTeX compilation
# Date: March 2018

# Data Wrangling
echo "Wrangling data ..."
python3 Wrangling.py
echo "Data wrangling complete"

# Model fitting
echo "Fitting models ..."
Rscript Model_fitting.R
echo "Model fitting complete"
echo "Plotting data ..."
Rscript Plotting_analysis.R
echo "Plotting complete"

rm *.pdf
# LaTeX
echo "Compiling LaTeX"
bash CompileLatex.sh ../Writeup/Miniproject.tex