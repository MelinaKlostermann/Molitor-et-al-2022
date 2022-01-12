#################################################################
#
#       RNA accessibility surronding PURA binding sites
#
#################################################################

### ===============================================================
### Structure of this script
### ---------------------------------------------------------------
###
### 1) Obtaining fasta files of bound sequences 
### 2) Obtaining fasta files of random background
### 3) Predicting RNA accessibility with RNAplfold (commandline) and concatinating output files
### 4) z-score calculation and plot
###
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

## PURA binding sites annotated on transcripts
BS_tx  # as obtaines from xxx


### ===============================================================
### 1) Obtaining fasta files of bound sequences 
### ===============================================================

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



### ================================================================
### 2) Obtaining fasta files of random background
### ================================================================

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




### ===============================================================
### 3) Predicting RNA accessibility with RNAplfold (commandline) and concatinating output files
### ================================================================

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




### ==========================================================
### 4) z-score calculation and plot
### ==========================================================


# 4.1) log-odds ratio
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

# 4.2) background mean and sd from subsets
#############################
### the population mean and background are calcualted from 1000 subsets of the background

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



# 4.3) z-scores of BS accesibility
#############################
# z-score  = sample mean -  population mean / poplulation sd (sample = binding site, population = background)

df_accessibility_bs_bg <- bg_stats_df_stats %>%
  mutate(bs_mean = gg_df_mean$BS, 
         z_score = (bs_mean - mean) /sd) %>% 
  rowwise() %>%
  mutate(p_z = 2*pnorm(-abs(z_score))) %>%
  ungroup() %>%
  mutate(p_z_adj = p.adjust(p_z),
         pos = -250:250,
         sig_p = p_z_adj < 0.05)



# 4.4) plot
#################
ggplot(gg_df_mean_z, aes(x=pos))+
  geom_hline(yintercept = 1, color = "grey")+
  geom_line(aes(y=z_score),  color = "black")+
  xlim(c(-100,100))+
  geom_rug(aes(color=sig_p))+
  scale_color_manual(values = c( "transparent", "darkgreen"))+
  theme(legend.position="none")





