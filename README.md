# substitutants_manuscript
Code for cancer datasets substitutant analysis

Cancer cells experience tryptophan shortage as a consequence of interferon-gamma (IFN) signaling pathway. It has been recently demonstrated that such tryptophan shortage leads to tryptophan to phenylalanine codon reassignment (W>F) resulting into aberrant peptide products. ABPEPserver offers the visualization of a large-scale proteomics analysis of multiple human cancer types in which aberrant peptides are detected. In this analysis, tryptophan to phenylalanine (W>F) codon reassignment were found to be a vastly abundant phenomenon in multiple cancer types. Furthermore, these W>F mis-incorporations, called W>F substitutants, were found to be enriched in tumors as compared with tumor-adjacent normal tissues, and often their appearance was associated with T-cell and oncogenic signaling activities. Proteomic cancer data from multiple cancer types, hosted by PDC commons database, have been used in this analysis. An in depth description of how the analysis is given in github link given in point 6.



Large Scale analysis of CPTAC exemplified by LSCC analysis

1. Datasets for LSCC are downloaded in MZML format from PDC commmons webpage https://proteomic.datacommons.cancer.gov/pdc/study/PDC000234

2. Datasets are formatted per TMT channel (eg, 01CPTAC) as instructed on the philosopher webpage https://github.com/Nesvilab/philosopher/wiki/Pipeline-mode-for-TMT-analysis

3. Philospher parameters are set according to Supplementary Tables in the manuscript.

3. Philospher was run;

bin/philosopher pipeline --config params/philosopher.yml 01CPTAC.....
(philosopher.yml parameters are submitted in supplementary tables, Database was created by substituting W-> to all other amino acids

4. ratio_peptide_None.tsv (output from philosopher) is used for further downstream analysis.

Extracted substitutants for all W>X combinations using UNIX command.
grep -Fwf W_F.txt ratio_peptide_None.tsv > W_F2.txtratio_peptide_None.tsv

5. Make a Rscript to get counts;

This script gives substitutant peptide filters, counts and plots (attached in the github folder).

peptide_filters.R (attached)


6. The detailed description in toto is here;

https://github.com/jasminesmn/ABPEPserver
