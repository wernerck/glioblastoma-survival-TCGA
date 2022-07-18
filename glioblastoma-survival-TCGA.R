# GBM survival analysis by low/high expression 
# UCSC Xena TCGA data
# Created in R 4.0.2

library(UCSCXenaTools)
library(dplyr)
library(survival)
library(survminer)
library(cgdsr)
library(magrittr)
library(ggplot2)
library(knitr)

# Download 3 tables using UCSCXenaTools: clinical, phenotype, rna
# Clinical table for OS and OS.time from TCGA (PANCAN)
cl = XenaGenerate(subset = XenaHostNames =="pancanAtlasHub")
cl1 = XenaFilter(cl, filterDatasets = "s1")

q1 = XenaQuery(cl1)
d1 = XenaDownload(q1)
clinical = XenaPrepare(d1)

# Phenotype, which has treatment info (original TCGA)
ph <- XenaGenerate(subset = XenaHostNames =="tcgaHub")
ph1 <- XenaFilter(ph, filterDatasets = "Clinical")
ph2 <- XenaFilter(ph1, filterDatasets = "TCGA.GBM.sample")

q2 <- XenaQuery(ph2)
d2 <- XenaDownload(q2)
phenotype <- XenaPrepare(d2)

# Rna, gene expression data for GBM patients (original TCGA microarray)
# Lets look at EGFR, a commonly studied gene
# 539 samples
rn <- XenaGenerate(subset = XenaHostNames =="tcgaHub")
rn1 <- XenaFilter(rn, filterDatasets = "GBM")
rn2 <- XenaFilter(rn1, filterDatasets = "ht_hg")

q3 <- XenaQuery(rn2)
d3 <- XenaDownload(q3)
rna <- XenaPrepare(d3)

# dim(rna) = 12042 (genes) x 539(samples) + 1(gene label)
# save the short version of the sample names
shortNames = names(rna)[2:540]
shortNamesF <- data.frame(shortNames) #converts list to data frame

# convert rna row to data frame; transpose to column
EGFRInd <- which(rna$sample == "EGFR")
EGFRData <- rna[EGFRInd,]
EGFRF <- data.frame(EGFRData[2:length(EGFRData)])
EGFRF <- t(EGFRF)
EGFR <- cbind(shortNames, EGFRF) #to make shortNames a column
colnames(EGFR) <- c("shortNames", "EGFR_exp") #to rename

# Data Cleaning
# filter clinical table by shortNames 
# leaves 518 patients
clin2 <- data.frame(clinical)
table1 <- merge(clin2, shortNamesF, by.x="sample", by.y="shortNames")

# Filter phenotype table by sample names
phen2 <- data.frame(phenotype)
table2 <- merge(phen2, table1, by.x="sampleID", by.y="sample")

# merge relevant patients into one large table: table1, table2, table3, EGFR
# yields matrix 518x197
data <- merge(table1, table2, by.x="sample", by.y="sampleID")
data2 <- merge(data, EGFR, by.x="sample", by.y="shortNames")

#create new table with only the columns we need 
new_data <- select(data2, sample, OS.x, OS.time.x, CDE_chemo_alk, CDE_radiation_standard, CDE_radiation_standard_probable, CDE_alk_chemoradiation_standard, CDE_therapy, EGFR_exp) %>% 
  rename(time = OS.time.x, status = OS.x, chemo = CDE_chemo_alk, rt = CDE_radiation_standard, rtp = CDE_radiation_standard_probable, chemort = CDE_alk_chemoradiation_standard, therapy = CDE_therapy)

new_data$EGFR_exp <- as.double(new_data$EGFR_exp)


# Define EGFR low/high cutoffs for gene expression

# zscore option
# new_data$EGFR_exp <- scale(new_data$EGFR_exp)
# new_data$prmt5_exp <- scale(new_data$prmt5_exp)

#define by quantile instead
EGFR_quant <- quantile(new_data$EGFR_exp, probs = c(0.25,0.5,0.75))

# Look at treatment type
# treatment: no treatment, chemotherapy, chemoradiation

# filter by treatments and prepare for analysis
nt_all <- filter(new_data, is.na(therapy)) %>% 
  select(sample, EGFR_exp, time, status)

chemo_all <- filter(new_data, chemo == "TRUE") %>% 
  select(sample, EGFR_exp, time, status)

chemort_all <- filter(new_data, rt == "TRUE" | rtp == "TRUE" | chemort == "TRUE")%>% 
  select(sample, EGFR_exp, time, status)


# Perform survival analysis by treatment and high/low expression of EGFR
categories <- c("nt_all", "chemo_all", "chemort_all")

for( val in seq(1, length(categories), 1)) {
  varName <- categories[val]
  
  varData0 <- eval(parse(text = varName))
  
  varDatam <- varData0 %>% 
    mutate(group = case_when(
      EGFR_exp > EGFR_quant[[3]] ~ 'EGFR High',
      between(EGFR_exp, EGFR_quant[[1]], EGFR_quant[[3]]) ~ 'EGFR Int',
      EGFR_exp < EGFR_quant[[1]] ~ 'EGFR Low',
      TRUE ~ NA_character_
    ))	
  
  fit_EGFR <- survfit(Surv(time, status) ~ group, data = varDatam)
  EGFR_plot <- ggsurvplot(fit_EGFR, font.x = 20, font.y = 20, font.title = c(24, "bold"), legend = c(0.8,0.8), title = paste(varName, "_EGFR", sep=""))
  plotName <- paste(varName, "_EGFR_plot.pdf", sep="")

  ggsave(plotName)}