## load packages
library(dplyr)

## load the yob data into a dataframe
id_yob  <- read.table("#redacted/data/raw/phenotypes/id_yob.txt", header = TRUE)

## load the SES variables
ses_df <- read.table("#redacted/data/raw/phenotypes/income_yrsedu.txt", header = TRUE) %>% rename(IID=eid)

## load the extra covariates
center_batch <- read.table("#redacted/data/raw/phenotypes/id_center_batch.txt", header = TRUE)

## load the panukb covariates
panukb_covs1 <- read.table("#redacted/data/raw/panukb/panukb_ancestry.txt", header = TRUE)
## remove extra covariates
panukb_covs <- panukb_covs1 %>% select(colnames(panukb_covs1)[1:11], sex)

## merge dfs
yob_covs1 <- merge(panukb_covs, id_yob, by = "IID")
yob_covs <- merge(yob_covs1, center_batch, by = "IID")
all_vars <- merge(yob_covs, ses_df, by = "IID")

## save
write.table(all_vars, "#redacted/data/raw/phenotypes/all_phenotypes.txt", row.names = FALSE, quote = FALSE)
