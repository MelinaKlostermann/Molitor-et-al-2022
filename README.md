# Molitor-et-al-2022

We analysed PURA's in vivo binding behaviour to RNAs in iCLIP experiments. Theses scripts contain the R code used to

* 1 Binding site definition and characterisation

+ Input: crosslinks in bw format, pureclip sites, Gencode Annotaiton (v31) as gtf
+ Output: crosslinks (bw_merges.rds), crossliks per sample (bw_all_samples.rds), binding sites (binding_sites.rds), filtered annotation (annotation.rds)

* 2 RNAseq changes in PURA knockdown
** Input: htseq count table, filtered annotation (annotation.rds)
** Output: deseq anaysis (deseq.rds), count table with tpms (rnaseq_counts.rds)



* Predict RNA accessibility and 5-mer frequencies considering mature transcript sequences
* Compare the crosslink patterns in three different PURA iCLIP experiments 


# Glossary

- crosslink site: number of nucleotide that has been found crosslinked to PURA once or several times
- crosslink events: number of crosslink events detected (one crosslink site can contain several crosslink events)
- binding site: PURA binding site (5nt wide) as defined in the binding site definition below

# See also

- Busch et al. 2019 - "iCLIP data analysis: A complete pipeline from sequencing reads to RBP binding sites" https://doi.org/10.1016/j.ymeth.2019.11.008
- Krakau, S., Richard, H. & Marsico, A. 2017 - "PureCLIP: capturing target-specific protein–RNA interaction footprints from single-nucleotide CLIP-seq data." https://doi.org/10.1186/s13059-017-1364-2
