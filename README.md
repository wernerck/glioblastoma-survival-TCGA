# Glioblastoma Survival by Low/High Gene Expression

TCGA (PANCAN) Gene expression data were downloaded from the University of California Santa Cruz cancer browser using the UCSCXenaTools package. Affymetrix U133a microarray data were selected for analysis, rather than RNAseq, because more microarray samples were available. Expression, treatment annotations (from the original TCGA database), and outcomes were available for 539 patients. 

Samples were stratified by treatment annotation into groups by treatment: no treatment, chemotherapy, and chemoradiation. Expression of the specified gene was split into "low" and "high" by median expression, and survival curves were plotted using the survminer package.

EGFR is used as the example gene in the code.