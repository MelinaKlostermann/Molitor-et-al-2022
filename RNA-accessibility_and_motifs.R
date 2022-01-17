#################################################################
#
#       RNA accessibility surronding PURA binding sites
#
#################################################################

### ===============================================================
### Structure of this script
### ---------------------------------------------------------------
###
### 1) Transfere binding sites to transcript annotation
### 2) Accessibility
##### 2.1) Obtaining fasta files of bound sequences 
##### 2.2) Obtaining fasta files of random background
##### 2.3) Predicting RNA accessibility with RNAplfold (commandline) and concatinating output files
##### 2.4) z-score calculation and plot
### 3) Motifes
### 3.1) 
###



### ===============================================================
### 0) Libraries & Input
### ===============================================================

# Libraries
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(tidyverse)
library(ggpubr)
library(Biostrings)

# Input
## transcript seqeunces fasta from gencode
transcript_annotation <- readDNAStringSet("./gencode.v31.transcripts.fa.gz")

## PURA binding sites 
BS_tx  # as obtaines from xxx

### ============================================================
### 1) map BS to transcripts
### ============================================================
###################
# get sequences of  mature transcripts
##################

# expressed transcripts ( = transcripts with any crosslinks)
annotation_transcripts <- annotation[annotation$type == "transcript"] 
expressed_transcripts <- subsetByOverlaps(annotation_transcripts, cl_endo) 

annotation_transcripts_exons <- annotation[annotation$type != "gene"] 
expressed_transcripts_GR_list <- subsetByOverlaps(annotation_transcripts_exons, expressed_transcripts) %>%
  splitAsList(., f = .$transcript_id) %>%
  GRangesList(.) 


# seq of all transcripts
transcript_annotation <- readDNAStringSet("/Users/melinaklostermann/Documents/projects/anno/GENCODEv31-p12/gencode.v31.transcripts.fa.gz") #alternative extractTranscriptSeqs(x, transcripts, ...)
transcript_anno_meta <- names(transcript_annotation) 
transcript_anno_meta <- data.frame(all = transcript_anno_meta) %>%
  tidyr::separate(., col = all,
                  into = c("transcript_id", "gene_id", "a", "b", "isoform_name", "gene_name", "entrez_gene_id", "gene_type"), sep = "\\|")


names(transcript_annotation) <- transcript_anno_meta$transcript_id


###########################
# BS seq considering mature transcripts
##########################


# prepare a txdb of expressed transcripts
expressed_transcripts_txdb <- unlist(expressed_transcripts_GR_list)
txdb <- makeTxDbFromGRanges(unlist(expressed_transcripts_GR_list))

# prepare a transcript mapper (contains transcript ids and names together with genomic positions of transcripts)
transcripts_txdb_mapper <- transcripts(txdb)

# get transcript realtive coordinates of BS
BS_tx<- mapToTranscripts(BS, txdb, extractor.fun = GenomicFeatures::exonsBy)


# change the seqnames to the trnascript names
BS_tx<- as.data.frame(BS_tx)
BS_tx$seqnames<- transcripts_txdb_mapper$tx_name[as.numeric(BS_tx$seqnames)] %>% substring(.,1,15)





### ============================================================
### 2) Accessibility
### ============================================================

#############################################
### 2.1) Obtaining fasta files of bound sequences 
#############################################

# get transcript_id and transcript lengths from fasta names
transcript_anno_meta <- names(transcript_annotation) 
transcript_anno_meta <- data.frame(all = transcript_anno_meta) %>%
  tidyr::separate(., col = all, into = c("transcript_id", "gene_id", "a", "b", "isoform_name", "gene_name", "entrez_gene_id", "gene_type"), sep = "\\|")
names(transcript_annotation) <- substring(transcript_anno_meta$transcript_id,1,15)
transcript_annotation_df <- data.frame(tx_name = names(transcript_annotation), width = width(transcript_annotation))


# elongate binding sites to window for RNAplfold analysis
w <- 248 

BS_tx_501nt <- BS_tx %>%
  left_join(transcript_annotation_df, by= c(seqnames = "tx_name"), suffix = c(".bs", ".tx")) %>% 
  mutate(end = end + w, start = start -w) %>%
  dplyr::filter((end <  width.tx)  & (start > 0))

BS_tx_501nt <- makeGRangesFromDataFrame(BS_tx_501nt, keep.extra.columns = T)
BS_tx_501nt <- BS_tx_501nt[width(BS_tx_501nt) == 501]

# only one transcript per BS
BS_tx_501nt <- BS_tx_501nt[!duplicated(BS_tx_501nt$xHits)]

# get transcript seqences of enlarged BS 
BS_tx_501nt_seqs <- getSeq(x = transcript_annotation, names = BS_tx_501nt)
BS_tx_501nt_seqs <- BS_tx_501nt_seqs[width(BS_tx_501nt_seqs)==501]

writeXStringSet(BS_tx_501nt_seqs, filepath = paste0(output,"endo_BS_trans_500nt.fasta"))



############################################
### 2.2) Obtaining fasta files of random background
############################################

set.seed(2)
n_random <- 10000

# filter for long enough transcripts
big_transcript_annotation_df <- transcript_annotation_df[transcript_annotation_df$width > 500,]

# get set of random transcripts
random_transcripts <- data.frame(transcript = sample(1:NROW(big_transcript_annotation_df), n_random, replace = F))
random_transcripts$transcript_id <- big_transcript_annotation_df[random_transcripts$transcript,]$tx_name

# seq of random expressed transcripts set
random_transcript_seqs <-  transcript_annotation[names(transcript_annotation) %in% random_transcripts$transcript_id,] 

# subset for random window of 300 nt
w_2 <- (2*w + 5)
random_transcript_seqs_500 <- random_transcript_seqs
random_transcript_seqs_pos <- list()

for(i in 1:NROW(random_transcript_seqs)){
  print(width(random_transcript_seqs[i])-(w_2-1))
  random_pos = sample(1:(width(random_transcript_seqs[i])-(w_2-1)), 1) 
  random_transcript_seqs_500[i] <- subseq( random_transcript_seqs[i], start = random_pos, width = w_2)
  
  random_transcript_seqs_pos[[i]] <- data.frame(idx = i, tx_name = names(random_transcript_seqs[i]), start = random_pos)
} 


writeXStringSet(random_transcript_seqs_500, filepath = paste0(output,"background_transcripts_500nt.fasta"))




###############################################
### 2.3) Predicting RNA accessibility with RNAplfold (commandline) and concatinating output files
###############################################

### Note RNAplfold is used here via command line on the fasta files created above
### RNAplfold creates on folder per fasta input
### each folder contains a text file ending on _*lunp for each sequence in the fasta file

input_folders <- list("./RNApl_background_transcripts_500nt/",
                      "./RNApl_endo_BS_trans_500nt/")

output_probs <- list("./RNApl_endo_BS_trans_500nt.RData")

# function to get probs
import_unpair_prob <- function(folder){
  RNApl_files =list.files(path = folder,
                          pattern = "*_lunp")
  
  unpaired_prob = lapply(RNApl_files, function(x) data.table::fread(paste0(folder,x)))
  
  names(unpaired_prob) = 1:length(unpaired_prob)
  
  unpaired_prob = unpaired_prob %>%
    map(~dplyr::rename(.x,  pos = `#i$`, `1`=`l=1`) %>%
          as.data.frame() ) %>%
    map_dfc(~.x[,2])
  
  return(unpaired_prob)
}

# make RData with probs
probs <- map(input_folders, ~import_unpair_prob(.x))
print("# probs")
head(probs)


save(probs, file =output_probs)




####################################
### 2.4) z-score calculation and plot
####################################


# 2.4.1) log-odds ratio
##############################
### First unpaired probabilitties from RNAplfold are tranfered to log-odd ratios 
### to obtain a bell-shaped distribution of values which is neccesary for z-score calculation


# log-odds ratio binding sites (BS)
accesibility_BS <- probs[[2]] %>% as.matrix(.)
accesibility_BS <- log(accesibility_BS/(1-accesibility_BS)) 
accesibility_BS[is.infinite(accesibility_BS)] <- NA
accesibility_BS <- as.data.frame(accesibility_BS)

# log-odds ratio background (bg)
accesibility_bg <- probs[[1]] %>% as.matrix(. )
accesibility_bg <- log(accesibility_bg/(1-accesibility_bg)) 
accesibility_bg[is.infinite(accesibility_bg)] <- NA
accesibility_bg <- as.data.frame(accesibility_bg)

# 2.4.2) background mean and sd from subsets
#############################
### the population mean and background are calculated from 1000 subsets of the background

# index for 1000 sets with 100 nt sequences
idx_sets <- replicate(1000, sample(1:ncol(accesibility_bg), 1000), simplify = F)

# get sets
bg_sets <- map(idx_sets, ~accesibility_bg[,.x])

# calc nt wise mean per set an make mean df
bg_stats <- map(bg_sets, ~t(.x) %>%
                  as.data.frame(.) %>%
                  summarise_all(.funs = function(x) mean(x, na.rm=T)))
bg_stats_df <- map_dfr(bg_stats, ~.x)

# calc mean and sd of means_df
bg_stats_df_stats <- data.frame(mean = apply(bg_stats_df, 2, function(x) mean(x, na.rm = T)),
                                sd = apply(bg_stats_df, 2, function(x) sd(x, na.rm =T)), pos = -250:250)



# 2.4.3) z-scores of BS accessibility
#############################
# z-score  = sample mean -  population mean / population sd (sample = binding site, population = background)

df_accessibility_bs_bg <- bg_stats_df_stats %>%
  mutate(bs_mean = gg_df_mean$BS, 
         z_score = (bs_mean - mean) /sd) %>% 
  rowwise() %>%
  mutate(p_z = 2*pnorm(-abs(z_score))) %>%
  ungroup() %>%
  mutate(p_z_adj = p.adjust(p_z),
         pos = -250:250,
         sig_p = p_z_adj < 0.05)



# 2.4.4) plot
#################
ggplot(gg_df_mean_z, aes(x=pos))+
  geom_hline(yintercept = 1, color = "grey")+
  geom_line(aes(y=z_score),  color = "black")+
  xlim(c(-100,100))+
  geom_rug(aes(color=sig_p))+
  scale_color_manual(values = c( "transparent", "darkgreen"))+
  theme(legend.position="none")







### ===========================
### Motif - 5-mer analysis
### ==========================

############################
# Sequences of 3'UTR and CDS BS
###########################
BS_3UTR <- BS_tx[BS_tx$BS_region == "3UTR"]
seqences_BS_3utr <- getSeq(x = transcriptome_seqs, names = BS_3UTR)

BS_CDS <- BS_tx[BS_tx$BS_region == "CDS"]
seqences_BS_cds <- getSeq(x = transcriptome_seqs, names = BS_CDS)


#############################
# random background 3'UTR
#############################
# get 3'UTRs on transcripts
threeUTR_transcript_anno <- mapToTranscripts(unlist(threeUTRsByTranscript(txdb, use.names = T)), txdb, extractor.fun = GenomicFeatures::exonsBy, use.names =T) %>%
  GenomicRanges::reduce() %>%
  as.data.frame() %>% 
  mutate(seqnames = substring(seqnames,1,15), strand ="+") %>%
  distinct(seqnames, .keep_all = T) %>%
  makeGRangesFromDataFrame() 

threeUTR_transcript_anno <- threeUTR_transcript_anno[seqnames(threeUTR_transcript_anno) %in% substring(expressed_transcripts$tx_name,1,15)] 

w_b <- 5

# select transcripts with crosslinks
r_transcripts = threeUTR_transcript_anno[sample(1:NROW(threeUTR_transcript_anno), replace = T)]
  
r_seqs = r_transcripts  %>%
  as.data.frame(.) %>%
  rowwise() %>% 
    mutate(start = sample(start:(end-w_b-2), 1),
           end = start+(w_b -1),
           width = 5) %>%
    ungroup() %>%
    as.data.frame() %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
  
r_seqs  = r_seqs[r_seqs %within% r_transcripts]

# random sampling  
random_bg_3utr = r_seqs[sample(1:NROW(r_seqs), NROW(seqences_BS_3utr))] 

# get sequence
seqences_bg_3utr <- getSeq(x = transcriptome_seqs, names =random_bg_3utr)


###################
# random background CDS
##################
CDS_transcript_anno <- mapToTranscripts(unlist(cdsBy(txdb, use.names = T)), txdb, extractor.fun = GenomicFeatures::exonsBy, use.names =T) %>%
  GenomicRanges::reduce() %>%
  as.data.frame() %>% 
  mutate(seqnames = substring(seqnames,1,15), strand ="+") %>%
  distinct(seqnames, .keep_all = T) %>%
  makeGRangesFromDataFrame() 

CDS_transcript_anno <- CDS_transcript_anno[seqnames(CDS_transcript_anno) %in% substring(expressed_transcripts$tx_name,1,15)] 

w_b <- 5


# select transcripts with crosslinks
r_transcripts = CDS_transcript_anno[sample(1:NROW(CDS_transcript_anno), replace = T)]
  
r_seqs = r_transcripts  %>% as.data.frame(.) %>% rowwise() %>% 
    mutate(start = sample(start:(end-w_b-2), 1),
           end = start+(w_b -1),
           width = 5) %>%
    ungroup() %>%
    as.data.frame() %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
# random sampling   
r_seqs  = r_seqs[r_seqs %within% r_transcripts]
  
# get sequence
random_bg_cds = r_seqs[sample(1:NROW(r_seqs), NROW(seqences_BS_cds))] 

# get sequences
seqences_bg_cds <- getSeq(x = transcriptome_seqs, names =random_bg_cds)


##########################
# frequency scatter
##########################
# function to calc oligo frequency
calc_olig_freq <- function(seq_regions, oligo){
  
  nu_freq = oligonucleotideFrequency(seq_regions, width = oligo) %>%
    colSums() %>% 
    as.data.frame()
  
  colnames(nu_freq) = c("nu_freq")
  nu_freq$olig_nu = rownames(nu_freq)
  nu_freq$rel_freq = nu_freq$nu_freq /sum(nu_freq$nu_freq)
  
  return(nu_freq)
}

# function to make 5 scatter plots of oligofrequency 
plot_freq_scatter <- function(gr_regions, bg, oligo, title, filename, window, color_thresh){
  
  # sites: get upstream and downstream window 
  gr_regions = gr_regions[(gr_regions+window+1) %within% transcriptome_seqs_coords]
  seq_regions = GRangesList(upstream = GRanges(seqnames = seqnames(gr_regions), ranges = IRanges(start = start(gr_regions) - window -1, width = window), strand = "+"),
                            BS = gr_regions,
                            downstream = GRanges(seqnames = seqnames(gr_regions), ranges = IRanges(start = end(gr_regions) + 1, width = window), strand = "+"))
  seq_regions = lapply(seq_regions, function(x) getSeq(transcriptome_seqs, names = x))
  
  # background: get upstream and downstream window 
  bg_regions = bg[(bg+window+1) %within% transcriptome_seqs_coords]
  seq_regions_bg = GRangesList(upstream = GRanges(seqnames = seqnames(bg_regions), ranges = IRanges(start = start(bg_regions) - window -1, width = window), strand = "+"),
                               BS = bg_regions,
                               downstream = GRanges(seqnames = seqnames(bg_regions), ranges = IRanges(start = end(bg_regions) + 1, width = window), strand = "+"))
  seq_regions_bg = lapply(seq_regions_bg, function(x) getSeq(transcriptome_seqs, names = x))
  
  
  # freqs
  freqs <- map(seq_regions, ~calc_olig_freq(.x, oligo = oligo))
  freqs_bg <- map(seq_regions_bg, ~calc_olig_freq(.x, oligo = oligo))  
  
  freqs <- map2(freqs, freqs_bg, ~left_join(.x, .y, by = "olig_nu", suffix = c(".bs", ".bg")))
  
  title_list = list("in binding site", paste(window, "nt upstream of binding site"), paste(window, "nt downstream of binding site"))
  freqs <- map2(freqs, title_list, ~ mutate(.x, group = .y,
                                            purines = stringr::str_count(olig_nu, pattern = c("G")) + stringr::str_count(olig_nu, pattern = c("A"))))
  freqs <- map_dfr(freqs, ~.x)
  freqs$group <- factor(freqs$group, levels = unlist(title_list)[c(2,1,3)])
  
  
  # make plots

  ggplot(freqs, aes(x = rel_freq.bg, y = rel_freq.bs, label = olig_nu, color = purines > color_thresh))+ 
    geom_point()+
    geom_point(data = freqs[freqs$purines > color_thresh,], aes(x = rel_freq.bg, y = rel_freq.bs, label = olig_nu), color = "blue")+
    theme_classic()+
    ggtitle(title)+
    geom_abline(color = "grey")+
    ggrepel::geom_text_repel(size = 2, max.overlaps = 10)+
    scale_color_manual(values = c("darkgrey", "blue"))+
    facet_wrap(~group, scales = "free_y")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(filename = filename,  width = 12, height=3, device = "pdf", path = outpath)
  
  
}

# use function to make plots
plot_freq_scatter(gr_regions = c(BS_CDS, BS_3UTR)+5, bg = c(random_bg_cds, random_bg_3utr), oligo = 5, window =20, color_thresh = 3,
                  filename= "all_BS+5_5nt_Freq_20nt_against_bg.pdf",
                  title="all_BS+5_5nt_Freq_20nt_against_bg")
