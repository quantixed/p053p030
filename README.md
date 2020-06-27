# p053p030
Code for p053p030

## Data

- `sequences` directory contains the DNA sequences of the plasmids used introducing GFP-FKBP tags during gene editing.
- `progression` has mitotic progression data for all cell lines



## Fiji Code

`knocksideways_spindle.ijm` will generate the data from a pre and post rapamycin image from a knocksideways experiment movie. Suggested data organisation: each movie has a directory with two images, run the code on this directory. These directories can be grouped by protein-of-interest and then by cell line (as described below.

## Igor Code

- `MitoticProgression.ipf` this code will read data from Excel workbooks. Data should be organised as follows: separate worksheets for the cells to be compared. Three columns per sheet labelled `NEB`, `Metaphase`, `Anaphase`. Each row is contains the frame number each cell reached these stages. Data used in the paper can be found in `Data/progression`

