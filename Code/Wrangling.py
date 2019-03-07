#!usr/bin/env python3
""" Data Wrangling for Miniproject """
__author__= 'Eva Linehan (eva.linehan18@imperial.ac.uk)'
__version__ = 0.01
__date__ = 'Jan 2019 '
__licence__ = 'CMEE Assessment'

import pandas as pd
import numpy as np

f = pd.read_csv('../Data/BioTraits.csv') # Import

f = f[f['OriginalTraitValue'] > 0] # Remove all Trait values that are less that 0
f = f[pd.notnull(f['OriginalTraitValue'])] # Remove NA's in all Trait values
f = f[pd.notnull(f['ConTemp'])] # Remove NA's in all Temperature measurements
points = 6
subset = f.groupby("FinalID").size() # Group ID's by number of measurements
subset = subset[subset >= points] # Create a subset of ID's with 6 observations or more
subset = subset.index.tolist() # Convert to list
f = f[f['FinalID'].isin(subset)] # Let new dataframe consist of FinalID's with those in the converted list
#f.shape

f['Kelvin'] = f['ConTemp'] + 273.15 # Create new column of Kelvin Temperature values

f.to_csv('../Data/Ready_to_fit.csv') # Export
