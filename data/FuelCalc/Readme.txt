This folder contains files to be used with FuelCalc software /https://www.firelab.org/project/fuelcalc-canopy-fuel-calculator-and-model). The goal is to produce plot-level canopy fuel estimates.

- "input_treelist.tre": tree lists
- "input_treelist.ffi": plot lists
These files are outputs from a previous analysis step. ALS point clouds are segmented into individual trees and reduced to tree-level metrics. File structure is modified so they can be used directly in FuelCalc to produce plot-level estimates.

- "species.csv"
This file contains species codes and calculation instructions for FuelCalc. It defines which allometric equations should be applied to calculate certain fuel components and needs to be placed in FuelCalc's root directory! 

- "test_file.ffi": template for creating plot list file.

