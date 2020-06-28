# p053p030
Code and data for Ryan et al. manuscript`#p053p030`

## Data

- `fixed_ksw` contains the outputs from Fiji (in `Fiji_outputs`) and the resulting data frames (in `R_combined`) to generate the plots in the paper.
- `progression` has mitotic progression data for all cell lines
- `sequences` directory contains the DNA sequences of the plasmids used introducing GFP-FKBP tags during gene editing.



## Fiji Code

`knocksideways_spindle.ijm` will generate the data from a pre and post rapamycin image from a knocksideways experiment movie. Suggested data organisation: each movie has a directory with two images, run the code on this directory. These directories can be grouped by protein-of-interest and then by cell line (as described below.

## Igor Code

- `SpindleQuantFromFiji.ipf` this code will read all the outputs from `knocksideways_spindle.ijm`. Data should be organised as follows: cell lines in separate directories, sub-directories for each protein assessed by live knocksideways. Consistent naming is required. Point Igor at the the directory which contains the cell line directories.
- `MitoticProgression.ipf` this code will read data from Excel workbooks. Data should be organised as follows: separate worksheets for the cells to be compared. Three columns per sheet labelled `NEB`, `Metaphase`, `Anaphase`. Each row is contains the frame number each cell reached these stages. Data used in the paper can be found in `Data/progression`


## R Code

- `ProcessFijiDataset.R` is a script to process the outputs from spindle measurements in Fiji. The images are blinded and require unblinding with a `log.txt` file in the same directory as the files.

- `fixed_ksw.R` is a script to process combined dataframes that are outputs from `fiji_ksw.R`. All data frames for one cell line  should be in `.rds` format and grouped in a directory named according to the cell line.

A further file, `lookup.csv` is required to run the code.