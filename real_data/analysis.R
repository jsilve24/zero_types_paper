library(zinbwave)
library(scRNAseq)
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(GenomicRanges)
library(stats4)
library(BiocGenerics)
library(parallel)
library(phyloseq)
library(tidyverse)
library(ggalt)
library(cowplot)
library(ggrepel)
library(ggridges)
library(broom)
library(MicrobeDS)
library(DESeq2)
library(parathyroidSE)
library(pasilla)
library(ALDEx2)
library(readxl)

set.seed(343421)

setwd("~/Research/zero_types/results/2019-07-11_real_data/")

# utility functions -------------------------------------------------------

load_pollen <- function(){
  data("fluidigm")
  
  # # Remove lowly expressed genes as in vignette
  # filter <- rowSums(assay(fluidigm)>3)>3
  # table(filter)
  # fluidigm <- fluidigm[filter,]
  
  # Take 100 most variables genes
  # assay(fluidigm) %>% log1p %>% rowVars -> vars
  # names(vars) <- rownames(fluidigm)
  # vars <- sort(vars, decreasing = TRUE)
  # head(vars)
  # fluidigm <- fluidigm[names(vars)[1:1000],]
  
  # just take biocondition 
  fluidigm <- fluidigm[,colData(fluidigm)$Biological_Condition %in% c("NPC" , "GW21")]
  Y <- assay(fluidigm)
  filter <- rowSums(Y>3)>3
  table(filter)
  Y <- Y[filter,]
  X <- model.matrix(~Biological_Condition + Coverage_Type, data=data.frame(fluidigm@colData))
  return(list(Y=Y, X=t(X), group_names=c("GW21", "NPC")))
}

load_zheng <- function(){
  # only 10X data seems to be properly formated by sonenson - could not load others
  # data.path <- "~/Research/data/_data_derived/sonenson/Soneson_datasets/GSE60749-GPL13112.rds"
  data.path <- "~/Research/data/_data_derived/sonenson/Soneson_datasets/10XMonoCytoT.rds"
  dat <- readRDS(data.path)
  dat <- updateObject(dat)
  
  # Using rounded length-scaled TPM as used in Sonenson
  Y <- assays(experiments(dat)[["gene"]])[["count_lstpm"]]
  pretty_gene_names <- function(g) paste0("G", substr(g, 10, 15))
  rownames(Y) <- pretty_gene_names(rownames(Y))
  
  # as in risso, do some basic filtering of really low abundance genes
  filter <- rowSums(Y>3)>3
  table(filter)
  Y <- Y[filter,]
  
  X <- cbind(1, dat$group=="cytotoxict")
  return(list(Y=Y, X=t(X), group_names=c("monocyte","cytotoxict")))
}

load_kostic <- function(){
  #Genomic analysis identifies association of Fusobacterium with colorectal carcinoma. 
  #Kostic, A. D., Gevers, D., Pedamallu, C. S., Michaud, M., Duke, F., Earl, A. M., et al. 
  #(2012). Genome research, 22(2), 292-298.

  # Following phyloseq vignette
  filepath = system.file("extdata", "study_1457_split_library_seqs_and_mapping.zip", package="phyloseq")
  kostic = phyloseq::microbio_me_qiime(filepath)
  kostic <- subset_samples(kostic, DIAGNOSIS != "None")
  kostic <- prune_samples(sample_sums(kostic) > 500, kostic)
  kostic <- filter_taxa(kostic, function(x) sum(x>3)>3, TRUE)
  
  
  Y <- as(otu_table(kostic), "matrix")
  rownames(Y) <- paste0("OTU",rownames(Y))
  
  # Already < 1000
  # Y %>% log1p %>% rowVars -> vars
  # names(vars) <- rownames(Y)
  # vars <- sort(vars, decreasing = TRUE)
  # Y <- Y[names(vars)[1:1000],]
  

  X <- rbind(1, as.numeric(sample_data(kostic)$DIAGNOSIS)-1)
  group_names <- c("Healthy", "Tumor")
  return(list(Y=Y, X=X, group_names=group_names))
}

load_gevers <- function(){
  data("RISK_CCFA")
  dat <- RISK_CCFA %>% 
    subset_samples(disease_stat!="missing", 
                   immunosup!="missing") %>% 
    subset_samples(diagnosis %in% c("no", "CD")) %>% 
    subset_samples(steroids=="false") %>% 
    subset_samples(antibiotics=="false") %>% 
    subset_samples(biologics=="false") %>% 
    subset_samples(biopsy_location=="Terminal ileum") %>% 
    prune_samples(sample_sums(.) >= 500,.) %>%
    filter_taxa(function(x) sum(x > 3) > 3, TRUE)
  Y <- as(otu_table(dat), "matrix")
  rownames(Y) <- paste0("OTU",rownames(Y))
  
  Y %>% log1p %>% rowVars -> vars
  names(vars) <- rownames(Y)
  vars <- sort(vars, decreasing = TRUE)
  Y <- Y[names(vars)[1:1000],]
  
  sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>% 
    mutate(age = as.numeric(as.character(age)),
           diagnosis = relevel(diagnosis, ref="no"), 
           disease_stat = relevel(disease_stat, ref="non-inflamed"))
  X <- t(model.matrix(~diagnosis, data=sample_dat))
  group_names <- c("Healthy", "CD")
  return(list(Y=Y, X=X, group_names=group_names))
}

# load_pasilla <- function(){ # Difficulty parsing this dataset
#   # as in DESeq2 vignette
#   pasCts <- system.file("extdata",
#                         "pasilla_gene_counts.tsv",
#                         package="pasilla", mustWork=TRUE)
#   pasAnno <- system.file("extdata",
#                          "pasilla_sample_annotation.csv",
#                          package="pasilla", mustWork=TRUE)
#   cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
#   coldata <- read.csv(pasAnno, row.names=1)
#   coldata <- coldata[,c("condition","type")]
#   rownames(coldata) <- sub("fb", "", rownames(coldata))
#   cts <- cts[, rownames(coldata)]
#   Y <- cts[rowSums(cts>3)>3,]
#   
#   # Y %>% log1p %>% rowVars -> vars
#   # names(vars) <- rownames(Y)
#   # vars <- sort(vars, decreasing = TRUE)
#   # Y <- Y[names(vars)[1:5000],]
# 
#   X <- model.matrix(~condition+type, data=coldata)
#   group_names <- c("treated", "untreated")
#   return(list(Y=Y, X = t(X), group_names=group_names))
# }

load_mcmurrough <- function(){
  data("selex")
  Y <- as.matrix(selex)
  X <-  rbind(1, c(rep(0, 7), rep(1, 7)))
  rownames(Y) <- unlist(map(strsplit(rownames(Y), ":"), function(x) paste(x, collapse="")))
  group_names <- c("NS", "S")
  return(list(Y=Y, X=X, group_names=group_names))
}

load_haglund <- function(){
  data("parathyroidGenesSE")
  Y <- assay(parathyroidGenesSE)
  coldata <- colData(parathyroidGenesSE)
  pretty_gene_names <- function(g) paste0("G", substr(g, 10, 15))
  rownames(coldata) <- coldata$run
  colnames(Y) <- coldata$run
  rownames(Y) <- pretty_gene_names(rownames(Y))
  filter <- rowSums(Y>3)>3
  Y <- Y[filter,]
  filter <- as.character(coldata$treatment) %in% c("Control", "DPN")
  filter <- colnames(Y)[filter]
  Y <- Y[,filter]
  coldata <- coldata[filter,]
  coldata$treatment <- as.character(coldata$treatment)
  X <- model.matrix(~treatment+time, data=coldata)
  group_names <- c("DPN", "Control")
  return(list(X=t(X), Y=Y, group_names=group_names))
}


fitZI_zinbwave <- function(data){
  fit <- zinbFit(Y=data$Y, X=t(data$X))
  DE <- getBeta_mu(fit)[2,]
  DiffZI <- getBeta_pi(fit)[2,] 
  return(list(DE=DE, DiffZI=DiffZI, Y=data$Y, X=data$X, group_names=data$group_names))
}

fitNZI_zinbwave <- function(data){
  fit <- zinbFit(Y=data$Y, X=t(data$X), O_pi=matrix(-1000000, ncol(data$Y), nrow(data$Y)))
  DE <- getBeta_mu(fit)[2,]
  return(list(DE=DE, Y=data$Y, X=data$X, group_names=data$group_names))
}

summarize_fits <- function(fitZI, fitNZI){
  results <- data.frame("ZI" = fitZI$DE, "NZI" = fitNZI$DE, 
                    "DiffZI" = fitZI$DiffZI, 
                    "gene" = rownames(fitZI$Y)) %>% 
    mutate(DiffDE = NZI-ZI)
  return(results)
}

identify_focus <- function(summary_fits){
  focus <- summary_fits %>% 
    arrange(-abs(DiffDE)) %>% 
    .[1:10,] %>% 
    .$gene %>% 
    as.character()
  summary_fits %>% 
    filter(gene%in%focus) %>% 
    arrange(-DiffDE) %>% 
    .$gene %>% 
    as.character()
}

plot_distributions <- function(data, focus){
  p <- data.frame(t(data$Y[focus,]), "Group" = data$group_names[data$X[2,]+1]) %>% 
    gather(gene, count, -Group) %>% 
    mutate(count = count+1, 
           gene = factor(gene, levels=rev(focus))) %>% 
    ggplot(aes(x = count, y=gene, fill=Group, color=Group)) +
    stat_binline(rel_min_height=0.001, alpha=0.7, draw_baseline = FALSE, scale=0.97, bins = 20) +
    scale_x_log10() +
    #theme_ridges() +
    theme_minimal()+
    xlab("Observed Counts") +
    ylab("Gene") +
    scale_fill_manual(values=c("#4daf4a", "#984ea3")) +
    scale_color_manual(values=c("#4daf4a", "#984ea3")) +
    #scale_fill_brewer(palette="Set1")+
    theme(legend.title = element_blank(), 
          axis.title.y = element_blank(), 
          text = element_text(size=13)) 
  return(p)
}

plot_rank <- function(fitZI, fitNZI){
  gene.names <- rownames(fitZI$Y)
  ZI.ranked <- gene.names[order(abs(fitZI$DE), decreasing=TRUE)]
  NZI.ranked <- gene.names[order(abs(fitNZI$DE), decreasing=TRUE)]
  ranks <- cbind(ZI.ranked, NZI.ranked)
  
  k.limit <- 50
  res <- matrix(0,k.limit, 2)
  colnames(res) <- c("Overlap", "K")
  for (k in 1:k.limit){
    res[k,2] <- k
    res[k,1] <- length(intersect(ranks[1:k,1], ranks[1:k, 2]))
  }
  p <- as.data.frame(res) %>% 
    ggplot(aes(x = K, y=Overlap)) +
    geom_point()+
    geom_path() +
    geom_segment(x = 0, y=0, xend=k.limit, yend=k.limit, color="red") +
    ylab("Genes in Common") +
    xlab("K") +
    ylim(c(0,k.limit))+
    theme_bw()
  return(p)
}

plot_DE <- function(summary_fits){
  focus <- identify_focus(summary_fits)
  p <- summary_fits %>% 
    mutate(label = ifelse(gene %in% focus, as.character(gene), NA )) 
  if (!is.null(summary_fits$DiffZI)) p <- mutate(p, `Differential\nZero-Inflation` = DiffZI)
  p <- ggplot(p, aes(x=NZI, y=ZI))
  if (!is.null(summary_fits$DiffZI)){
    p <- p + geom_jitter(aes(color=`Differential\nZero-Inflation`), alpha=1)
  } else {
    p <- p + geom_jitter(alpha=1)
  }
  p <- p + 
    theme_minimal()+
    xlab("Log Differential Expression (NB)") +
    ylab("Log Differential Expression (ZINB)") +
    geom_text_repel(aes(label=label), size=2.5)
  if (!is.null(summary_fits$DiffZI)){
   p <- p+ scale_color_gradient2(mid="grey", high="#e41a1c", low="#377eb8", midpoint = median(p$data$DiffZI))
  }
  return(p)
}

plot_combined <- function(summary_fits, data, focus){
  p1 <- plot_DE(summary_fits)
  p2 <- plot_distributions(data, focus)
  plot_grid(p1, p2, align = "h")
}

# analysis ----------------------------------------------------------------

run_and_save <- function(data, fn){
  fitZI <- fitZI_zinbwave(data)
  fitNZI <- fitNZI_zinbwave(data)
  summary_fits <- summarize_fits(fitZI, fitNZI)
  focus <- identify_focus(summary_fits)
  
  # Create rank plot
  p <- plot_rank(fitZI, fitNZI)
  ffn <- paste0("rank_",fn,".pdf")
  ggsave(ffn, plot=p, height=3, width=3, units="in")
  
  # Create combined plot
  p <- plot_combined(summary_fits, data, focus)
  ffn <- paste0("combined_", fn, ".pdf")
  ggsave(ffn, plot=p, height=4, width=8, units="in")
  
  # A few key statistics
  if (!is.null(summary_fits$DiffZI)){
    test <- broom::tidy(cor.test(summary_fits$DiffZI, summary_fits$DiffDE, method="spearman"))
    write.table(test, file=paste0("DiffZI_DiffDE_corr_",fn,".tsv"))
  }
}

summarize_datasets <- function(l){
  ns <- c("Sparsity")
  summary_fxn <- function(x){
    out <- rep(0, 1)
    names(out) <- ns
    out["Sparsity"] <- sum(x==0)/prod(dim(x))
    return(out)
  }
  
  out <- l %>% 
    map("Y") %>% 
    map(summary_fxn) %>% 
    do.call(rbind,.)
  colnames(out) <- ns
  write.table(out, file="datasets_summary.tsv")
}


# apply code --------------------------------------------------------------

run_and_save(load_pollen(), "pollen")
run_and_save(load_zheng(), "zheng")
run_and_save(load_kostic(), "kostic")
run_and_save(load_gevers(), "gevers")
run_and_save(load_mcmurrough(), "mcmurrough")
run_and_save(load_haglund(), "haglund")

l <- list()
l$pollen <- load_pollen()
l$zheng <- load_zheng()
l$kostic <- load_kostic()
l$gevers <- load_gevers()
l$mcmurrough <- load_mcmurrough()
l$haglund <- load_haglund()
summarize_datasets(l)
