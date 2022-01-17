#################################################################
#
#       Binding site definition
#
#################################################################

### ===============================================================
### Structure of this script
### ---------------------------------------------------------------
###
### 1) Preprocess input files
### 2) Definition of binding sites
### 3) Characterization of binding sites
### 4) Transfer of binding sites to transcript coordinates
###
###
###



### ===============================================================
### 0) Libraries & Input
### ===============================================================
library(GenomicRanges)
library(rtracklayer)
library(knitr)
library(GenomicFeatures)
library(dplyr)
library(ggpubr)
library(Gviz)
library(biomaRt)

#######################
# Input
########################
# bigwig files of crosslink events (all 4 samples merged)
bw_all_plus_path 
bw_all_minus_path 

# bigwig files of crosslink events (all 4 samples separate)
bw_1_plus_path 
bw_1_minus_path 

bw_2_plus_path 
bw_2_minus_path 

bw_3_plus_path 
bw_3_minus_path 

bw_4_plus_path 
bw_4_minus_path 


# pureclip calls (obtained by running pureclip on pseudo samples 1u2 and 3u4 see below)
pureclip_path  

# gene annotation (gencode annotation v31)
mygft 


### ===============================================================
### 1) Preprocess input files
### ===============================================================

################################
# Split crosslinks to 1nt events
################################
split_bw_crosslinks_to_1_nt <- function(bw){
  # split ranges in 1 nt events
  bw_split <- exomeCopy::subdivideGRanges(bw, subsize=1)
  # match scores by an overlap index
  idx <- findOverlaps(bw_split, bw)
  # add scores
  bw_split$score <- bw[subjectHits(idx)]$score
  # readd strand info
  strand(bw_split) <- strand(bw)
  return(bw_split)
}

bw_all_samples <- list(bw_1_plus, bw_1_minus, bw_2_plus, bw_2_minus, bw_3_plus, bw_3_minus, bw_4_plus, bw_4_minus) %>% lapply(., function(x) split_bw_crosslinks_to_1_nt(bw=x))

bw_merges <- list(bw_all_plus, bw_all_minus) %>% lapply(., function(x) split_bw_crosslinks_to_1_nt(bw=x))

########################################
# Load and filter Genecode annotation
#######################################

filter_gft_anno <- function(gft_anno, standard_chrom, protein_coding, include_GL_3, longest_transcript){
  library(GenomicRanges)
  library(dplyr)
  library(GenomicFeatures)
  
  #load gft
  gft <- rtracklayer::import(gft_anno)
  
  # Standard Chromosomes
  if(standard_chrom == T){
    gft <-keepStandardChromosomes(gft, pruning.mode = "coarse")
  }
  
  # Protein Coding Transcripts
  if(protein_coding == T){
    gft <- gft[gft$gene_type=="protein_coding"]
    gft <- c(gft[gft$type!="gene" & gft$transcript_type=="protein_coding"], gft[gft$type=="gene"])
  }
  
  # Gene level 1 or 2
  gft_GL<- gft[gft$level <= 2]
  
  # Gene level 3
  if(include_GL_3==T){
    gft_GL3 <- gft[gft$level==3 & !(gft$gene_id %in% gft_GL$gene_id)]    
    gft_GL <- c(gft_GL, gft_GL3)  
  }
  
  # Transcriptlevel <=3 or NA
  gft_GL_TL <-gft_GL[!is.na(gft_GL$transcript_support_level) & gft_GL$transcript_support_level <= 3]
  gft_TL_NA <- gft_GL[is.na(gft_GL$transcript_support_level)]
  gft_Transcripts <- c(gft_GL_TL, gft_TL_NA[!(gft_TL_NA$gene_id %in% gft_GL_TL$gene_id)])
  gft_gen_trans <- c(gft_Transcripts, gft_GL[gft_GL$type=="gene"]) # readd genes (have transcript support level NA)
  
  # longest transcript
  if(longest_transcript==T){
    # use genomicfeatures file to get information about transcript length (sum of exons)
    gft_GenFeat<- makeTxDbFromGRanges(gft_gen_trans)
    transcriptLength <- transcriptLengths(gft_GenFeat)
    
    # add transcript_length to granges
    gft_GL_TL_dframe <- as.data.frame(gft_gen_trans)
    gft_with_tx_len <- merge(transcriptLength[, c("tx_name", "tx_len")], gft_GL_TL_dframe, by.x ="tx_name", by.y = "transcript_id",
                             all.y = TRUE)
    
    gft_longestTranscript <- gft_with_tx_len[gft_with_tx_len$type == "transcript", ] %>%
      dplyr::group_by(.$gene_id) %>%
      dplyr::arrange(dplyr::desc(.$tx_len)) %>%
      dplyr::slice(1) %>%
      ungroup()
    
    gft_regions_longestTranscript <- gft_gen_trans[gft_gen_trans$transcript_id %in% gft_longestTranscript$tx_name]
    gft_region_genes_LT <- c(gft_regions_longestTranscript, gft_gen_trans[gft_gen_trans$type=="gene" & gft_gen_trans$gene_id %in% gft_regions_longestTranscript$gene_id])
    gft_gen_trans <- gft_region_genes_LT
  }
  return(gft_gen_trans)
}

annotation <- filter_gft_anno(gft_anno = mygft,
                              standard_chrom = T, 
                              protein_coding = F, 
                              include_GL_3 = T, 
                              longest_transcript = F)


########################
# get pureclip output
#########################
pureclip_sites <- import(pureclip_path, format = "bedgraph")

#clean up columns
pureclip_sites <- as.data.frame(pureclip_sites) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
pureclip_sites$NA.2 <- NULL
pureclip_sites$score <- pureclip_sites$NA.
pureclip_sites$NA. <- NULL
strand(pureclip_sites) <- pureclip_sites$NA.1
pureclip_sites$NA.1 <- NULL
pureclip_sites$round_score <- round(pureclip_sites$score, digits = 1)

pureclip_sites <- keepStandardChromosomes(pureclip_sites, pruning.mode = "coarse")
pureclip_sites


### =====================================
### 2) Definition of binding sites
### Note: Adapted from Busch et al. 2019 - "iCLIP data analysis: A complete pipeline from sequencing reads to RBP binding sites"
### ======================================

############################
# make 5nt binding sites
###########################
Define_Binding_Sites <- function(pureclip, bw_plus, bw_minus, windowsize, out){
  
  # Merge Gaps < 8 from single pureclip sites
  pureclip = GenomicRanges::reduce(pureclip, min.gapwidth = 8)
  
  #remove sites with 1 or 2 nt length
  pureclip = pureclip[width(pureclip) > 2]
  
  
  bw_plus = import.bw(bw_plus, as="Rle")
  bw_minus = import.bw(bw_minus, as= "Rle")
  
  
  
  final.peaks.plus.gr <- GRanges()
  final.peaks.minus.gr <- GRanges()
  
  
  #Initialize the remaining PureCLIP CL regions to check for peaks
  remaining.regions.plus.gr <- subset(pureclip, strand == "+")
  remaining.regions.minus.gr <- subset(pureclip, strand == "-")
  
  window.radius <- (windowsize-1)/2
  while(TRUE){
    
    #No regions left to check for peaks
    if (length(remaining.regions.plus.gr) == 0 & length(remaining.regions.minus.gr) == 0){
      break
    }
    
    if (length(remaining.regions.plus.gr) != 0 ){
      #Get the raw CL counts in the remaining PureCLIP CL regions
      # returns rle list of all regions and turns it into matrix
      raw.remaining.PureCLIP.CL.regions.plus.m <- as.matrix(bw_plus[remaining.regions.plus.gr])
      
      #Identify the center of the PureCLIP CL regions (Position with max counts)
      # and store its indice
      raw.remaining.PureCLIP.CL.regions.plus.m[
        is.na(raw.remaining.PureCLIP.CL.regions.plus.m)] <- -Inf # set Na to -infinite
      max.pos.indice.plus <- max.col(raw.remaining.PureCLIP.CL.regions.plus.m, 
                                     ties.method = "first")
      
      
      #Create a peak region of 9nt that is centered to the max position
      peaks.plus.gr <- remaining.regions.plus.gr
      start(peaks.plus.gr) <- start(peaks.plus.gr) + max.pos.indice.plus - 1
      end(peaks.plus.gr) <- start(peaks.plus.gr)
      peaks.plus.gr <- peaks.plus.gr + window.radius
      
      
      #Store the new peaks
      final.peaks.plus.gr <- c(final.peaks.plus.gr, peaks.plus.gr)
      
      #Remove the peaks from the CL regions to search for additional peaks
      #Excise additionally 4 nucleotides up and downstream
      peaks.plus.grl <- as(peaks.plus.gr+window.radius, "GRangesList")
      
      remaining.regions.plus.gr <- unlist(psetdiff(remaining.regions.plus.gr, peaks.plus.grl))
    }
    if (length(remaining.regions.minus.gr) != 0 ){
      #Get the raw CL counts in the remaining PureCLIP CL regions
      # returns rle list of all regions and turns it into matrix
      raw.remaining.PureCLIP.CL.regions.minus.m <- as.matrix(
        bw_minus[remaining.regions.minus.gr])
      
      #Identify the center of the PureCLIP CL regions (Position with max counts) 
      #and store its indice
      raw.remaining.PureCLIP.CL.regions.minus.m[
        is.na(raw.remaining.PureCLIP.CL.regions.minus.m)] <- -Inf
      max.pos.indice.minus <- max.col(raw.remaining.PureCLIP.CL.regions.minus.m, ties.method = "last")
      
      #Create a peak region of 9nt that is centered to the max position
      peaks.minus.gr <- remaining.regions.minus.gr
      start(peaks.minus.gr) <- start(peaks.minus.gr) + max.pos.indice.minus - 1
      end(peaks.minus.gr) <- start(peaks.minus.gr)
      peaks.minus.gr <- peaks.minus.gr + window.radius
      
      #Store the new peaks
      final.peaks.minus.gr <- c(final.peaks.minus.gr, peaks.minus.gr)
      
      #Remove the peaks from the CL regions to search for additional peaks
      #Excise additionally 4 nucleotides up and downstream
      peaks.minus.grl <- as(peaks.minus.gr+window.radius, "GRangesList")
      
      remaining.regions.minus.gr <- unlist(psetdiff(remaining.regions.minus.gr,
                                                    peaks.minus.grl))
    }
  }
  export(final.peaks.plus.gr, 
         con= paste(out,"_plus.bed", sep = ""), 
         format = "bed")
  export(final.peaks.minus.gr, 
         con= paste(out,"_minus.bed", sep = ""), 
         format = "bed")
  save(final.peaks.minus.gr, file= paste(out,"_plus.RData", sep = ""))
  save(final.peaks.plus.gr, file= paste(out,"_minus.RData", sep = ""))
  
  returnlist <- list(peaks.minus = final.peaks.minus.gr, peaks.plus = final.peaks.plus.gr)
  return(returnlist)
}


binding_sites <- Define_Binding_Sites(pureclip = pureclip_sites, 
                                      bw_plus = bw_all_plus_path, 
                                      bw_minus = bw_all_minus_path,
                                      windowsize = 5, # windowsize - size that Bindingsites should have 
                                      out =  "./Binding_site_windows_5nt" )



############################
# Keep only BS with PureCLIP center
############################
# get all BS
binding_sites <- c(final.peaks.minus.gr, final.peaks.plus.gr)
#get centers
BS_centers <- binding_sites - 2

#keep only overlaps with pureclip sites
pureclip_sites_anno_bw_SOB <- makeGRangesFromDataFrame(pureclip_sites, 
                                                       keep.extra.columns = TRUE)
binding_sites_center_PS <- binding_sites[queryHits(findOverlaps(
  BS_centers, pureclip_sites))]



###########################
# Keep only BS with max PureCLIP site at center
##########################
# get bw rles
bw_plus_rle <- import.bw(bw_all_plus_path, as="Rle")
bw_minus_rle <- import.bw(bw_all_minus_path, as="Rle")

# split BS by strand
binding_sites_center_PS_plus <- binding_sites_center_PS[strand(binding_sites_center_PS)=="+"]
binding_sites_center_PS_minus <- binding_sites_center_PS[strand(binding_sites_center_PS)=="-"]

# make matrix of BS
binding_sites_center_PS_plus_m <- as.matrix(bw_plus_rle[binding_sites_center_PS_plus])
binding_sites_center_PS_minus_m <- as.matrix(bw_minus_rle[binding_sites_center_PS_minus])

# calc max for each BS (one BS is one row in the matrix)
max_BS_plus <- apply(binding_sites_center_PS_plus_m,1,max)
max_BS_minus <- apply(binding_sites_center_PS_minus_m,1,max)

# subset for center = max
binding_sites_center_PSmax_plus <- binding_sites_center_PS_plus[
  max_BS_plus == binding_sites_center_PS_plus_m[,3]] 
binding_sites_center_PSmax_minus <- binding_sites_center_PS_minus[
  max_BS_minus == binding_sites_center_PS_minus_m[,3]]



###########################
# Keep only BS with at least 2  crosslink sites
############################
binding_sites_center_PSmax_plus_m <- as.matrix(bw_plus_rle[binding_sites_center_PSmax_plus])
binding_sites_center_PSmax_minus_m <- as.matrix(bw_minus_rle[binding_sites_center_PSmax_minus])

crosslink_sites_plus <- apply(binding_sites_center_PSmax_plus_m, 1, function(x) 5-sum(x == 0)) 
crosslink_sites_minus <- apply(binding_sites_center_PSmax_minus_m, 1, function(x) 5-sum(x == 0)) 


binding_sites_center_PSmax_plus_2cl <- binding_sites_center_PSmax_plus[crosslink_sites_plus > 1]
binding_sites_center_PSmax_minus_2cl <- binding_sites_center_PSmax_minus[crosslink_sites_minus > 1]




### ===================================
### 3) Characterization of binding sites
### ====================================

#######################################
### Sample Reproducibility 
#######################################
#get bws as rles
sample1.minus.rle <- import.bw( bw_1_minus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")
sample2.minus.rle <- import.bw( bw_2_minus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")
sample3.minus.rle <- import.bw( bw_3_minus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")
sample4.minus.rle <- import.bw( bw_4_minus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")

sample1.plus.rle <- import.bw( bw_1_plus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")
sample2.plus.rle <- import.bw( bw_2_plus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")
sample3.plus.rle <- import.bw( bw_3_plus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")
sample4.plus.rle <- import.bw( bw_4_plus_path, as="Rle") %>% keepStandardChromosomes(pruning.mode = "coarse")

# Sum up cl events per binding site 
load("/Users/melinaklostermann/Documents/projects/PURA/02_R_new_pip/01-BS_def/02-BS_def_endo/PURA_endo_R/Report-10-plots/Binding_site_windows_5nt_final.RData")
bs.p = binding_sites_final[strand(binding_sites_final) == "+"]
bs.p$clp_rep1 = sample1.plus.rle[bs.p] %>% sum
bs.p$clp_rep2 = sample2.plus.rle[bs.p] %>% sum
bs.p$clp_rep3 = sample3.plus.rle[bs.p] %>% sum
bs.p$clp_rep4 = sample4.plus.rle[bs.p] %>% sum


bs.m = binding_sites_final[strand(binding_sites_final) == "-"]
bs.m$clp_rep1 = sample1.minus.rle[bs.m] %>% sum
bs.m$clp_rep2 = sample2.minus.rle[bs.m] %>% sum
bs.m$clp_rep3 = sample3.minus.rle[bs.m] %>% sum
bs.m$clp_rep4 = sample4.minus.rle[bs.m] %>% sum

# Combine
binding_sites = c(bs.p, bs.m)

repro_df <- data.frame(s1 = binding_sites$clp_rep1,
                      s2 = binding_sites$clp_rep2,
                      s3 = binding_sites$clp_rep3,
                      s4 = binding_sites$clp_rep4) 

# plot
repro_scatter_df <- ecdf_df[,1:4] %>% mutate(s1 = log2(s1), s2 = log2(s2), s3=log2(s3), s4=log2(s4)) %>%
  mutate(s1 = case_when(s1== -Inf ~ 0, T ~ s1),
         s2 = case_when(s2== -Inf ~ 0, T ~ s2),
         s3 = case_when(s3== -Inf ~ 0, T ~ s3),
         s4 = case_when(s4== -Inf ~ 0, T ~ s4))

scatter_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    ggrastr::rasterise( geom_point(), dpi = 300) +
    ggrastr::rasterise(stat_density2d(aes(fill=..level..), geom="polygon"), dpi = 300) +
    scale_fill_gradientn(colours=report_color) +
    coord_cartesian(xlim = c(0,12.5), ylim = c(0,12.5))+
    geom_abline(slope=1, colour = "darkgrey", linetype="dashed")
  p
}

cor_fun <- function(data, mapping, method="pearson", ndp=2, sz=3, stars=T, ...){
  
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  corr <- cor.test(x, y, method=method)
  est <- corr$estimate
  lb.size <- sz* abs(est) 
  
  if(stars){
    stars <- c("***", "**", "*", "")[findInterval(corr$p.value, c(0, 0.001, 0.01, 0.05, 1))]
    lbl <- paste0(method, ": ", round(est, ndp), stars)
  }else{
    lbl <- round(est, ndp)
  }
  
  ggplot(data=data, mapping=mapping) + 
    annotate("text", label=lbl, x= 6, y= 6, size=lb.size,...)+
    theme(panel.grid = element_blank())
}


GGally::ggpairs(repro_scatter_df, upper = list(continuous = cor_fun), lower = list(continuous = scatter_fn), title = "Reproducability matrix - comparisons of 2 samples")+
  theme_bw()


#############################
# bound genes
#############################
BS_endo_genes <- annotation_genes_dup_rem[substr(annotation_genes_dup_rem$gene_id,1,15) %in% substr(BS_endo$gene_id,1,15)]  %>% as.data.frame()  %>% filter(!grepl(gene_type, pattern = "pseudo") & gene_type != "TEC" ) 

ggplot(BS_endo_genes, aes(x = gene_type))+
  geom_bar(fill = "blue")+
  stat_count(aes(label=paste0(sprintf("%1.11f", ..count../sum(..count..)*100),
                              "% \n", ..count..), y=0.5*max(..count..)), 
             geom="text", colour="black", size=4, position=position_dodge(width=1)) +
  coord_flip()+
  theme_paper()



################################################
# assign BS to their regions with a hierarchical approach
###############################################

#split up 3' and 5' UTR
# make a df with needed info from UTRs
utr_df <- as.data.frame(annotation[annotation$type=="UTR"])
#get cds and match start of cds to utr df
cds4UTR <- annotation[annotation$type=="CDS"]%>% as.data.frame
idx_UTR_cds <- match(utr_df$transcript_id, cds4UTR$transcript_id)
utr_df <- cbind(utr_df, start_CDS=cds4UTR[idx_UTR_cds,]$start)

# UTR is 3' if the start of the cds is downstream of the start of the UTR, else is 5'
utr_df <- utr_df %>% 
  mutate(type=ifelse((.$strand == '+'& .$start > .$start_CDS)|
                       (.$strand=='-' & .$start <.$start_CDS), "3UTR", "5UTR"))
utr_35 <- utr_df %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
utr_35$start_CDS <- NULL

annotation <- c(annotation[annotation$type!="UTR"], utr_35)


#make a GRanges with all regions overlaping BS
idx_BS_regions <- findOverlaps(binding_sites_repro, annotation, ignore.strand = FALSE)
BS_with_all_regions <- binding_sites_repro[queryHits(idx_BS_regions)]
elementMetadata(BS_with_all_regions) <- c(elementMetadata(BS_with_all_regions),
                                          elementMetadata(annotation[subjectHits(idx_BS_regions)]))

table(BS_with_all_regions$gene_type)

# hirarcy: 3'UTR, 5'UTR, CDS, intron 
# start and stop codon left out -> part of UTRs

# initialise column of peak region
BS_with_all_regions$BS_region <- NA

#########
# 3' UTR
#########
# get BS overlapping with 3'UTRs 
BS_3utr <- BS_with_all_regions[BS_with_all_regions$type == "3UTR"] %>% unique 
BS_3utr$BS_region <- "3UTR"
BS_assigned_regions <- BS_3utr


##########
# 5'UTR
#########
# get cds BS, that are not overlaping to utrs
BS_5utr <- BS_with_all_regions[-queryHits(findOverlaps(
  BS_with_all_regions, BS_assigned_regions, type="any"))] %>% .[.$type=="5UTR"]%>% unique
BS_5utr$BS_region <- "5UTR"
# add assigned CDS BS to assigned regions
BS_assigned_regions <- c(BS_assigned_regions, BS_5utr)

##########
# CDS
########
# get cds BS, that are not overlaping to utrs
BS_cds <- BS_with_all_regions[-queryHits(findOverlaps(
  BS_with_all_regions, BS_assigned_regions, type="any"))] %>% .[.$type=="CDS"]%>% unique
BS_cds$BS_region <- "CDS"
# add assigned CDS BS to assigned regions
BS_assigned_regions <- c(BS_assigned_regions, BS_cds)


###############
# non coding 
#############
# get non-coding BS, that are not overlaping to before
BS_noncod <- BS_with_all_regions[-queryHits(findOverlaps(
  BS_with_all_regions, BS_assigned_regions, type="any"))] %>% .[.$type=="exon"]%>% unique
BS_noncod$BS_region <- "non_cod"
# add assigned CDS BS to assigned regions
BS_assigned_regions <- c(BS_assigned_regions, BS_noncod)


########## 
# intron
##########
# get intron BS, that are not overlaping to before
BS_intron <- BS_with_all_regions[-queryHits(findOverlaps(
  BS_with_all_regions, BS_assigned_regions, type="any"))] %>% .[.$type=="transcript"]%>% unique
BS_intron$BS_region <- "intron"
BS_assigned_regions <- c(BS_assigned_regions, BS_intron)

#####################
# small RNAs are annotaed as genes
#####################
BS_other <- BS_with_all_regions[-queryHits(findOverlaps(
  BS_with_all_regions, BS_assigned_regions, type="any"))] %>% unique
BS_other$BS_region <- "non_cod"
BS_assigned_regions <- c(BS_assigned_regions, BS_other)

table(BS_assigned_regions$BS_region)






