# p053p030
Code and data for Ryan & Shelford et al. manuscript`#p053p030`

**Defining endogenous TACC3–chTOG–clathrin–GTSE1 interactions at the mitotic spindle using induced relocalization**

Ellis L Ryan^, James Shelford^, Teresa Massam-Wu, Richard Bayliss, Stephen J Royle

*bioRxiv* [2020.07.01.181818](https://doi.org/10.1101/2020.07.01.181818)

## Data

-  `live_ksw` contains the outputs from Fiji (`knocksideways_spindle.ijm`). Organised by cell line, these csv files can be used to generate ternary diagrams and arrow plots for live knocksideways experiments using `SpindleQuantFromFiji.ipf`.
- `progression` has mitotic progression data for all cell lines
- `sequences` directory contains the DNA sequences of the plasmids used introducing GFP-FKBP tags during gene editing.

Data for R plots are in the `R` directory.


## Fiji Code

`knocksideways_spindle.ijm` will generate the data from a pre and post rapamycin image from a knocksideways experiment movie. Suggested data organisation: each movie has a directory with two images, run the code on this directory. These directories can be grouped by protein-of-interest and then by cell line (as described below.

## Igor Code

- `SpindleQuantFromFiji.ipf` this code will read all the outputs from `knocksideways_spindle.ijm`. Data should be organised as follows: cell lines in separate directories, sub-directories for each protein assessed by live knocksideways. Consistent naming is required. Point Igor at the the directory which contains the cell line directories.
- `MitoticProgression.ipf` this code will read data from Excel workbooks. Data should be organised as follows: separate worksheets for the cells to be compared. Three columns per sheet labelled `NEB`, `Metaphase`, `Anaphase`. Each row is contains the frame number each cell reached these stages. Data used in the paper can be found in `Data/progression`


## R Code

Three R projects to generate plots in the paper.

- `fixed_ksw_figure` 
- `gtse1_lidl_figure`
- `gtse1_fragment_figure`

All work in a similar way. Data (outputs from Fiji) are in the `Data` directory and can be processed experiment-by-experiment or cell line by cell line to generate dataframes that are saved to `Output/Data`. A further script then calculates the statistics and plots the data.