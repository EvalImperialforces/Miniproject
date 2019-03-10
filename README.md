# Title: Miniproject
# Author: Eva Linehan el1718@ic.ac.uk
# Date: March 2019
# Course : MSc Computational Methods in Ecology and Evolution


## Execute Miniproject

### Clone repository
To execute this Miniproject, the repository must first be copied to a local machine and saved in a chosen directory. This can be downloaded or cloned using 'git clone' in a terminal followed by the SSH key.

### Prerequisites
This project requires the following software versions and related packages in order to be executable;

1. Python version 3.6:

	+ NumPy (1.13.3)
	+ pandas (0.23.4)

2. R version 3.4.4:

	+ dplyr (0.8.0.1)
  	+ ggplot2 (3.1.0)
  	+ gridExtra (2.3)
  	+ Hmisc (4.2.0)
  	+ minpack.lm (1.2.1)
  	+ plyr (1.8.4)
  	+ RColorBrewer (1.1.2)
  	+ reshape2 (1.4.3)
  	+ stringr (1.4.0)
  	+ xtable (1.8.3)

3. Bash 4.4.19

4. pdfTeX version 3.14 


### Running the project
From the command line, navigate to the '**Code**' directory within this Miniproject directory and use bash to execute the 'run_Miniproject.sh' script.
This will wrangle the BioTraits.csv dataset in the '**Data**' folder, using 'Wrangling.py', to produce 'Ready_to_fit.csv' used for model fitting. The model fitting script, 'Model_fitting.R' will be timed and produce 'Model_Results.csv' for plotting. All plots and figures, generated from 'Plotting_analysis.R' will be stored in the '**Results**' directory. 
Finally the script will compile 'Miniproject.tex' and 'References.bib' to produce a pdf version of the final report in the '**Writeup**' directory.




## Folder contents
This folder contains CMEE coursework relating to the Miniproject submission. All project work is arranged in the '**Data**', '**Code**', '**Results**' and '**Writeup**' folders.

The '**Data**' folder comprises of data files used for analysis. '**Code**' contains code used to execute data wrangling, analysis and plotting. '**Results**' contains all plots and figures generated as part of the Miniproject assessment and '**Writeup**' contains .tex files and related materials used in producing the final report.


### '**Code**' contains the following; 

1. 'run_MiniProject.sh'file
This shell script runs all code and LaTeX scripts to reproduce the Miniproject from start to finish. Bash version 4.4.19 is used in running this script.

2. 'Wrangling.py'file
This python script imports 'BioTraits.csv' to filter and subset data for model fitting before exporting the data as 'Ready_to_fit.csv'. 
Python version 3.6.7 is used, including 'NumPy' and 'pandas' libraries for data manipulation. 

3. 'Model_fitting.R'file 
This R script imports 'Ready_to_fit.csv' generated from data wrangling, and fits models to the data. 
The 'compiler' library allows functions to be pre-compiled.
The package 'minpack.lm' is required to perform non-linear least squares fitting for non-linear models. Model outputs are exported as 'Model_Results.csv' for plotting and analysis. 
R version 3.4.4 is used.

4. 'Plotting_analysis.R'file
This R script imports 'BioTraits.csv' and 'Model_Results.csv' to plot data and results. 
The 'compiler' library allows functions to be pre-compiled.
The package 'minpack.lm' is required to perform non-linear least squares fitting when fitting models to plotted data.
The package 'ggplot2' is required to create plots along with packages 'RColorBrewer', 'stringr' and 'Hmisc' which allow for colour and label customization in LaTeX-readable syntax. 
The package 'gridExtra' and 'xtable' are used to customize and convert tables to .tex files for LaTeX import. 
Packages 'dplyr', 'plyr' and 'reshape2' are required to manipulate, aggregate and melt dataframes for data visualization.
All figures and tables are exported as .pdf or .tex files to the 'Results' folder.
R version 3.4.4 is used.

5. 'CompileLatex.sh' file
This bash script compiles .tex documents from LaTeX into .pdf format while performing a clean-up of assitional unwanted outputs.

### '**Data**' contains the following;

1. 'BioTraits.csv' file
Database containing standardized trait values, temperatures and other biologically relevant data that is wrangled before analysis.

2. 'Ready_to_fit.csv' file
The BioTraits dataset filtered and condensed, containing only data included in the analysis.

3. 'Model_Results.csv' file
Comma seperated file containing sample ID and related model outputs for plotting and analysis.


### 'Results' contains the following;

1. Models_fitted' folder
Contains 1,577 figures generated from 'Plotting_analysis.R' in which models are fitted to each individual sample.

2. 'Appendix2.pdf' file
.pdf file that contains figure 2 in the Appendix of the final report.

3. 'Fig_2.pdf'
.pdf file that contains figure 2 printed in final report.

4. 'Fig_3.pdf'
.pdf file that contains figure 3 printed in final report.

5. 'Fig_4.pdf'
.pdf file that contains figure 4 printed in final report.

6. 'Fig_5.pdf'
.pdf file that contains figure 5 printed in final report.

7. 'Samples_boxplot.pdf'
.pdf file that contains figure 1 in the Appendix of the final report.

8. 'Table1.tex'
.tex file that contains table 1 printed in final report.

9. 'Table2.tex'
.tex file that contains table 2 printed in final report.

10. 'Table3.tex'
.tex file that contains table 3 printed in final report.


### '**Writeup**' contains the following;

1. 'crest.png' file
.png image used on coverpage of final report.

2. 'Miniproject.tex' file
.tex file containing report write up to be read and compiled by LaTeX.

3. 'References.bib' file
.bib file containing bibliography for final report.

4. 'Miniproject.pdf' file
.pdf file containing compiled final report.
