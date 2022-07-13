## load packages
library(dplyr)
library(metafor)

## load the files into a dataframe
afr_roh  <- read.table("#redacted/data/processed/roh_estimates/afr_roh.hom.indiv", header = TRUE) %>% select(IID, NSEG, KB, KBAVG)
afr_funi  <- read.table("#redacted/data/processed/roh_estimates/afr_funi.ibc", header = TRUE) %>% select(IID, Fhat3)
amr_roh  <- read.table("#redacted/data/processed/roh_estimates/amr_roh.hom.indiv", header = TRUE) %>% select(IID, NSEG, KB, KBAVG)
amr_funi  <- read.table("#redacted/data/processed/roh_estimates/amr_funi.ibc", header = TRUE) %>% select(IID, Fhat3)
csa_roh  <- read.table("#redacted/data/processed/roh_estimates/csa_roh.hom.indiv", header = TRUE) %>% select(IID, NSEG, KB, KBAVG)
csa_funi  <- read.table("#redacted/data/processed/roh_estimates/csa_funi.ibc", header = TRUE) %>% select(IID, Fhat3)
eas_roh  <- read.table("#redacted/data/processed/roh_estimates/eas_roh.hom.indiv", header = TRUE) %>% select(IID, NSEG, KB, KBAVG)
eas_funi  <- read.table("#redacted/data/processed/roh_estimates/eas_funi.ibc", header = TRUE) %>% select(IID, Fhat3)
eur_roh  <- read.table("#redacted/data/processed/roh_estimates/eur_roh.hom.indiv", header = TRUE) %>% select(IID, NSEG, KB, KBAVG)
eur_funi  <- read.table("#redacted/data/processed/roh_estimates/eur_funi.ibc", header = TRUE) %>% select(IID, Fhat3)
mid_roh  <- read.table("#redacted/data/processed/roh_estimates/mid_roh.hom.indiv", header = TRUE) %>% select(IID, NSEG, KB, KBAVG)
mid_funi  <- read.table("#redacted/data/processed/roh_estimates/mid_funi.ibc", header = TRUE) %>% select(IID, Fhat3)

## load the phenotypic data into a dataframe
all_vars  <- read.table("#redacted/data/raw/phenotypes/all_phenotypes.txt", header = TRUE)
all_vars$center <- as.factor(all_vars$center)
all_vars$batch <- as.factor(all_vars$batch)

## run associations
# Vector of ancestries in data set
ancestries <- c("afr", "amr", "csa", "eas", "eur", "mid")
anc_name <- c("African", "Admixed American", "Central/South Asian", "East Asian", "European", "Middle Eastern")

# Run association analyses and table results
results_table1 <- NULL # initialize empty df
for (j in 1:length(ancestries)){
## load roh data
roh_df <- get(paste(ancestries[j],"_roh",sep=""))

## calc froh
roh_df$froh <- roh_df$KB/(2.77*10^6)

## load fgrm data
fgrm_df <- get(paste(ancestries[j],"_funi",sep=""))

## merge all fstats data
df1 <- merge(roh_df, fgrm_df, by = "IID")

## merge all data
df2 <- merge(df1, all_vars, by = "IID")
df2$sex <- as.factor(df2$sex)

## run main regressions
froh_main <- lm(scale(froh) ~ scale(yob) + sex + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) + scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + center + batch, data = df2)
funi_main <- lm(scale(Fhat3) ~ scale(yob) + sex + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) + scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + center + batch, data = df2)
## save results
results1 <- cbind("sample" = "UKB",
                  "ancestry" = paste(ancestries[j]),
                  "estimate" = (c(summary(froh_main)$coefficients[2,1])),
                  "se" = (c(summary(froh_main)$coefficients[2,2])),
                  "p" = (c(summary(froh_main)$coefficients[2,4])),
                  "model" = "froh_main",
                  "N"= (c(nobs(froh_main))))

results2 <- cbind("sample" = "UKB",
                  "ancestry" = paste(ancestries[j]),
                  "estimate" = (c(summary(funi_main)$coefficients[2,1])),
                  "se" = (c(summary(funi_main)$coefficients[2,2])),
                  "p" = (c(summary(funi_main)$coefficients[2,4])),
                  "model" = "funi_main",
                  "N"= (c(nobs(funi_main))))

## run regressions with ses covariates
froh_ses <- lm(scale(froh) ~ scale(yob) + scale(years_edu) + scale(income) + sex + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) + scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + center + batch, data = df2)
funi_ses <- lm(scale(Fhat3) ~ scale(yob) + scale(years_edu) + scale(income) + sex + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) + scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10) + center + batch, data = df2)
## save results
results3 <- cbind("sample" = "UKB",
                  "ancestry" = paste(ancestries[j]),
                  "estimate" = (c(summary(froh_ses)$coefficients[2,1])),
                  "se" = (c(summary(froh_ses)$coefficients[2,2])),
                  "p" = (c(summary(froh_ses)$coefficients[2,4])),
                  "model" = "froh_ses",
                  "N"= (c(nobs(froh_ses))))

results4 <- cbind("sample" = "UKB",
                  "ancestry" = paste(ancestries[j]),
                  "estimate" = (c(summary(froh_ses)$coefficients[2,1])),
                  "se" = (c(summary(froh_ses)$coefficients[2,2])),
                  "p" = (c(summary(froh_ses)$coefficients[2,4])),
                  "model" = "froh_ses",
                  "N"= (c(nobs(froh_ses))))


results_table1 <- as.data.frame(rbind(results_table1, results1, results2, results3, results4))
}


results_table2 <- results_table1


## display table
results_table2

## save table
write.table(results_table2, "#redacted/results/ukb_all_models.txt", row.names = FALSE, quote = FALSE)
