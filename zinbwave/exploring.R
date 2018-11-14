library(scRNAseq)
library(zinbwave)
library(tidyverse)
library(xtable)
library(ggridges)

setwd("~/Research/zero_types/results/2018-04-13_zinbwave/")


# Preprocess and Run Model ------------------------------------------------

data("fluidigm")

# Remove lowly expressed genes as in vignette
filter <- rowSums(assay(fluidigm)>5)>5
table(filter)
fluidigm <- fluidigm[filter,]

# Take 100 most variables genes
assay(fluidigm) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(fluidigm)
vars <- sort(vars, decreasing = TRUE)
head(vars)
fluidigm <- fluidigm[names(vars)[1:100],]

# just take biocondition 
fluidigm <- fluidigm[,colData(fluidigm)$Biological_Condition %in% c("NPC" , "GW21")]
x <- colData(fluidigm)$Biological_Condition

#methods(class=class(fitZI))
fitZI <- zinbFit(fluidigm, X="~Biological_Condition + Coverage_Type")
fitNZI <- zinbFit(fluidigm, X="~Biological_Condition + Coverage_Type",
                  epsilon_min_logit=100000000000000)


# Investigate Results -----------------------------------------------------

# Ensure low Pi in NZI
getLogitPi(fitZI) %>% max()
getLogitPi(fitNZI) %>% max()

### Look for Relationship between zirm and diff.DE ###
focus <- order(abs(getBeta_mu(fitZI)[2,]-getBeta_mu(fitNZI)[2,]), decreasing=TRUE)[1:10]
focus.signdiff <- which(sign(getBeta_mu(fitZI)[2,])*sign(getBeta_mu(fitNZI)[2,])==-1)
focus <- c(focus, focus.signdiff)
focus.gene.names <- colnames(t(assay(fluidigm)[focus,]))

ln.NPC.to.GW21.ZI <- getBeta_mu(fitZI)[2,]
ln.NPC.to.GW21.NZI <- getBeta_mu(fitNZI)[2,]
ln.NPC.to.GW21 <- data.frame("ZI"=ln.NPC.to.GW21.ZI, "NZI" = ln.NPC.to.GW21.NZI, 
                             "gene" = colnames(t(assay(fluidigm)))) %>% 
  gather(Model, label, -gene) %>% 
  mutate(label = signif(log2(exp(label)), 2)) %>% 
  spread(Model, label)

tmp <- t(assay(fluidigm)) %>%
  as.data.frame() %>%
  bind_cols(., "Biological_Condition" = colData(fluidigm)$Biological_Condition) %>%
  gather(gene, count, -Biological_Condition) %>%
  group_by(Biological_Condition, gene) %>%
  summarise("nZero" = sum(count==0),
            "n" = n(), 
            "NonZeroMean" = mean(count[count>0]),
            "Mean" = mean(count)) %>%
  arrange(gene) %>%
  ungroup() 
tmpG <- filter(tmp, Biological_Condition=="GW21") %>% 
  select(-Biological_Condition)
tmpN <- filter(tmp, Biological_Condition=="NPC") %>% 
  select(-Biological_Condition)

tmp <- full_join(tmpG, tmpN, by="gene", suffix=c(".GW21", ".NPC")) %>%
  mutate(zirm =  log2(Mean.NPC/Mean.GW21)-log2(NonZeroMean.NPC/NonZeroMean.GW21)) %>% 
  full_join(ln.NPC.to.GW21, by="gene") %>% 
  mutate(diff.DE = NZI-ZI)

tmp %>% 
  select(gene, zirm, diff.DE) %>% 
  arrange(zirm) %>%
  mutate(focus1 = gene %in% colnames(t(assay(fluidigm)[focus[1:10],])), 
         focus2 = gene %in% colnames(t(assay(fluidigm)[focus.signdiff,]))) %>% 
  mutate(focus = ifelse(focus1, "Max Diff", FALSE)) %>% 
  mutate(focus = ifelse(focus2, "Sign Diff", focus)) %>% 
  mutate(focus = ifelse(focus==FALSE, "Other", focus)) %>% 
  mutate(focus = factor(focus, levels=c("Max Diff", "Sign Diff", "Other"))) %>% 
  rename(Filter = focus) %>% 
  ggplot(aes(x=zirm, y = diff.DE)) +
  geom_point(aes(color=Filter)) + 
  theme_minimal() +
  ylab("Difference in Differential Expression\n(NB-ZINB)") +
  xlab("Zero Inclusion Ratio of Means") +
  scale_color_manual(values=c("red", "blue", "black"))
ggsave("ZIRM_visualization.pdf", height=5, width=8, units="in")

# test for relationship (among all transcripts)
cor.test(tmp$zirm, tmp$diff.DE, method="spearman")
length(tmp$zirm)


### MAKE TABLE ###

tab <- tmp %>% 
  mutate(nZero.NPC = paste0(nZero.NPC,"/", n.NPC), 
         nZero.GW21 = paste0(nZero.GW21,"/", n.GW21)) %>% 
  mutate("nZero\n(NPC,GW21)" = paste0("(", nZero.NPC, ",", nZero.GW21, ")")) %>% 
  mutate(rankNZI = rank(-NZI, ties.method="max"), 
         rankZI = rank(-ZI, ties.method="max")) %>% 
  rename(ZIRM = zirm) %>% 
  mutate(NZI = round(NZI, 2), 
         ZI = round(ZI, 2), 
         ZIRM = round(ZIRM, 2), 
         diff.DE = round(diff.DE, 2)) %>%
  mutate(ZI = paste0(ZI, " (", rankZI, ")"), 
         NZI = paste0(NZI, " (", rankNZI, ")")) %>% 
  select(gene, `nZero\n(NPC,GW21)`, NZI, ZI, diff.DE, ZIRM) %>% 
  as.data.frame() %>% 
  .[match(focus.gene.names, .$gene),]
rownames(tab) <- tab$gene
tab %>% 
  select(-gene) %>% 
  xtable()
  


### Make Plot ###
data.frame(t(assay(fluidigm)[focus,]), "Biological_Condition"=x) %>% 
  arrange(Biological_Condition) %>% 
  gather(gene, count, -Biological_Condition) %>% 
  mutate(count = count+1, 
         gene = factor(gene, levels=rev(focus.gene.names))) %>%
  ggplot(aes(x = count, y = gene, fill=Biological_Condition)) +
  stat_density_ridges(rel_min_height=0.01, alpha=0.7) +
  scale_x_log10() +
  theme_ridges() +
  xlab("Count") +
  ylab("Gene") +
  scale_fill_brewer(palette="Set1")
ggsave("diff_genes.pdf", height=7, width=8, units="in")
