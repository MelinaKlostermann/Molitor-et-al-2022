# Molitor-et-al-2022

We analysed PURA's in vivo binding behaviour to RNAs in iCLIP experiments. Theses scripts contain the R code used to

## 1 Binding site definition and characterisation (HeLa, endogenous PURA)

+ Input: crosslinks in bw format, pureclip sites, Gencode Annotaiton (v31) as gtf
+ Output: crosslinks (bw_merges.rds), crossliks per sample (bw_all_samples.rds), binding sites (binding_sites.rds), filtered annotation (annotation.rds)

## 2 RNAseq changes in PURA knockdown
* Input: htseq count table, filtered annotation (annotation.rds)
* Output: deseq anaysis (deseq.rds), count table with tpms (rnaseq_counts.rds)

## 3 Characterisation of PURA binding (HeLa, endogenous PURA)
* Input:  binding sites (binding_sites.rds), filtered annotation (annotation.rds), crosslinks (bw_merges.rds), gff3 of annotation, count table with tpms (rnaseq_counts.rds)
* Output: binding sites with additional infos (binding_sites_characterized.rds)

## 4 Comparison of the crosslink patterns in three different PURA iCLIP experiments from HeLa cells and NPC PURA iCLIP 
* Input: binding sites (binding_sites.rds), filtered annotation (annotation.rds), crosslink bws of all 4 iCLIP experiments

## 5 Predict RNA accessibility and 5-mer frequencies considering mature transcript sequences (HeLa, endogenous PURA)
* Input: fasta of transcripts,  binding sites with additional infos (binding_sites_characterized.rds), filtered annotation (annotation.rds), crosslinks (bw_merges.rds)

## 6 Comparison of binding to RNAseq changes
* Input: deseq anaysis (deseq.rds), binding sites with additional infos (binding_sites_characterized.rds), crosslinks (bw_merges.rds)

## 7 Shotgun proteomics
* Input: maxquant table, deseq anaysis (deseq.rds), count table with tpms (rnaseq_counts.rds),  binding sites with additional infos (binding_sites_characterized.rds)

## 8 Reviewer figure: difference in proteomics results between ProteinDiscoverer and Maxquant
* Input: maxquant table with matching runs, maxquant table without matching runs, proteomeDiscoverer table (thermofisher), peptide table maxquant (with matching), peptide table thermofisher

## 9 Overlap of PURA binding with dendritic transcriptome  and with to p-body and stress granule enriched RNAs
* Input: dendritic RNA table (Middleton et al 2019), RNAseq of pbody vs whole cell (Hubstenberger et al), RNAseq of stress granule vs whole cell (Khong et al),
deseq anaysis (deseq.rds), count table with tpms (rnaseq_counts.rds),  binding sites with additional infos (binding_sites_characterized.rds)







# Glossary

- crosslink site: number of nucleotide that has been found crosslinked to PURA once or several times
- crosslink events: number of crosslink events detected (one crosslink site can contain several crosslink events)
- binding site: PURA binding site (5nt wide) as defined in the binding site definition below

# See also

- Busch et al. 2019 - "iCLIP data analysis: A complete pipeline from sequencing reads to RBP binding sites" https://doi.org/10.1016/j.ymeth.2019.11.008
- Krakau, S., Richard, H. & Marsico, A. 2017 - "PureCLIP: capturing target-specific protein–RNA interaction footprints from single-nucleotide CLIP-seq data." https://doi.org/10.1186/s13059-017-1364-2
- RNApLfold: https://www.tbi.univie.ac.at/RNA/RNAplfold.1.html
- Zhou Y, Zarnack K (2022). cliProfiler: A package for the CLIP data visualization. R package version 1.2.0, https://github.com/Codezy99/cliProfiler. 
- Middleton et. al 2019 (doi:10.1186/s12915-019-0630-z)

