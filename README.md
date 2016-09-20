# BaP_microbialInteractions
Data and analysis of microbial interactions that change in the presence of benzopyrene.

## Images
The images of colonies at 2x are stored here with the naming schema <media_condition>_<strain1>-<strain2>_<replicate_ID>

media_condition can be either "cntrl" or "benzo".

strain can be "p14" (P. aeruginosa PA14), "p1" (P. aeruginosa PA01), "s" (S. aureus), "hi" (H. influenzae"), or "hp" (H. parainfluenzae)

replicate_ID ranges from 0001 to 0006.

## Analysis
The comma-delimited file "EditedObjects.csv" is the output from [CellProfiler](http://cellprofiler.org/) of the analysis of the images. 
The R script "analyzeColonySizes.R" contains all the code to analyze the data in "EditedObjects.csv".

## Figures
Contains figures produced by "analyzeColonySizes.R".

## Questions and Contact
Any questions can be directed to Matt Biggs (mb3ad [at] virginia [dot] edu).

