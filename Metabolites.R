
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(dplyr)
library(robustbase)
library(StatMatch)
library(sva) # BiocManager::install("sva")
library(RColorBrewer)
library(colourvalues)
library(readr)
library(Rcpm) # devtools::install_github("ricoderks/Rcpm") # Git needs to be installed
library(bnlearn)
library(Rgraphviz) # BiocManager::install("Rgraphviz")

# Import dataset ----------------------------------------------------------

df_metabolome <- read_csv("../rawdata/metabodata_npa_metab_id.csv") %>% arrange(uid) 


# Batch effect corection --------------------------------------------------

df_mapping <- read_csv("../rawdata/mapping_metab_id_name.csv") 

covadjdata <- readstata13::read.dta13("../rawdata/m35_metadata_modified2019.dta") 
batchdata <- read.csv(file="../rawdata/batch_npa.csv", header=TRUE, sep=",",check.names=FALSE)
phenodata <-  batchdata %>% merge(covadjdata,by="study_id",all.y=F)
phenodata1 <- phenodata[,-2]
rownames(phenodata1) <- phenodata[,2]

normdata <- read.csv(file="../rawdata/metabodata_npa_metab_id_forcombat_n1013.csv", row.names=1, header=TRUE)
normdatamat <- as.matrix(normdata)
batch = phenodata1$batch
modcombat = model.matrix(~IntensiveTreatment+RVonly, data=phenodata1) # "mod" include 1 or 2 variables associated with batch number.
combat_normdata1 = ComBat(dat=normdatamat, mod=modcombat, batch=batch, par.prior=TRUE, prior.plots=FALSE) # mod=modcombat
df_metabolome_combat<-t(combat_normdata1)  # data includes negative values. This limits available normalization methods.


# Data normalization ------------------------------------------------------

# PQN for searching top metabolites

df_metabolome_norm <- pqn(df_metabolome_combat)

# MinMax scaling for the Deepinsight analysis
MinMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

for (i in 1:285){
  df_metabolome_norm[,i] <- MinMax(df_metabolome_norm[,i])
}

# extract uid names -------------------------------------------------------

uid <- rownames(df_metabolome_norm)
df_metabolome_norm_uid <- data.frame(uid,df_metabolome_norm)
rownames(df_metabolome_norm_uid) <- NULL

# Import metadata and outcome (asthma) data -------------------------------

df_metadata <- readstata13::read.dta13("../rawdata/m35_metadata_modified2019.dta") %>%
  dplyr::select(uid, study_id, 
         Age_mo,intake_sex,raceethn,premature37,intake_child_weight_born,mata_delivery,
         prev_breathingprob,icu_previous,intake_daycare,intake_smoke,mata_smoke,
         parent_asthma, parent_eczema,intch_weight,intch_respiratory_rate,o2sat_initial,
         eos_4perc,IgE_any_food,CPAPintubate,IntensiveTreatment,LOS_dys,abx_pre_inp,
         corticosteroids_pre_inp,RSV,HRV,RVA_NPAseq,RVB_NPAseq,RVC_NPAseq,
         recurrent_wheeze_sleep_36mo,RWsleep_time2event_36mo,inpatient_hfo,
         birth_season,aeroallergens_atopy,PAM4,site,antibiotics_life,corticosteroids_life,
         CPAPintubate,IntensiveTreatment) %>% tibble::as_tibble() %>%
  mutate(severity = ifelse(CPAPintubate == 1, 1, ifelse(inpatient_hfo == 1, 1, 0)))

# metadata6yasthma <- read.csv("M35_Asthma_6yr_2021_0223_forKH.csv", header=TRUE)
# df_metadata = merge(metadata6yasthma,df_metadata,by="study_id",all.y=T)


# Re-coding metadata ------------------------------------------------------

df_metadata$raceethn <- df_metadata$raceethn %>% as.factor()
df_metadata$birth_season <- df_metadata$birth_season %>% as.factor()
df_metadata$prev_breathingprob2 <- ifelse(df_metadata$prev_breathingprob==0, 0, 1) %>% as.factor()
df_metadata$parent_asthma[is.na(df_metadata$parent_asthma)] <- 0
df_metadata$parent_eczema[is.na(df_metadata$parent_eczema)] <- 0
df_metadata$anyIgE <- ifelse(df_metadata$IgE_any_food==1, 1,
                             ifelse(df_metadata$aeroallergens_atopy==1, 1,
                                    0)) %>% as.factor()
df_metadata$RVA_NPAseq[is.na(df_metadata$RVA_NPAseq)] <- 0
df_metadata$RVB_NPAseq[is.na(df_metadata$RVB_NPAseq)] <- 0
df_metadata$RVC_NPAseq[is.na(df_metadata$RVC_NPAseq)] <- 0


# Merge metabolites data and metadata -------------------------------------

metabo_variables <- merge(df_metabolome_norm_uid,df_metadata,by="uid") # metabolome data
metab_data <- metabo_variables[,c(2:286,288:290,325)]

# MWAS --------------------------------------------------------------------

variable <- NULL
for (i in 1:285){ 
  severe.fit <- glm(severity ~ ., 
                    family = binomial(),
                    metab_data[,c(i,286:289)]) # covariates:intake_sex, Age_mo, raceethn
  sum <- summary(severe.fit)
  sum_met <- sum$coefficients[2,] %>% as.data.frame() %>% t()
  rownames(sum_met) <- colnames(metab_data)[i]
  variable <- rbind(variable,sum_met)
}

var_select <- filter(as.data.frame(variable), variable[,4]<0.05)

# Merge metabolites and transcriptome data --------------------------------

transcriptome_variables <- read_csv("../transcriptome.csv")
metab_data <- metabo_variables %>%  dplyr::select("study_id",rownames(var_select)) 
fwrite(metab_data,"../metabolites.csv")

metab_trans <- merge(metab_data, transcriptome_variables, by = "study_id")
metab_trans <- metab_trans[,c(1,979:984,2:978)]

fwrite(metab_trans,"../trans_metab.csv")

# Make interaction term ---------------------------------------------------

for (i in 8:135){ # This loop takes about 15 minutes
  for(k in 136:984){
    add_column <- data.frame(metab_trans[,i] * metab_trans[,k])
    colnames(add_column) <- paste0(colnames(metab_trans)[i],"-",colnames(metab_trans)[k])
    metab_trans <- cbind(metab_trans, add_column)
  }
} 

## MinMax scaling for interaction terms

for (i in 985:109656){
  metab_trans[,i] <- MinMax(metab_trans[,i])
}

metab_trans_int <- metab_trans
fwrite(metab_trans_int,"../trans_metab_int.csv")

# (After the main analysis) Bayesian network analysis ---------------------
## https://qiita.com/hrkz_szk/items/a213c2c4ba823cbf78f6

feature_severe <- read_csv("../feature_severe.csv")
bayesian_var <- met_tra_variables[,c("severity",colnames(feature_severe))]

dag = hc(bayesian_var, score = "bic-g") # whitelist = wl, blacklist = bl
dag

graphviz.plot(dag, shape = "ellipse")

## boot strap

set.seed(1)
str.diff = boot.strength(bayesian_var, R = 200, algorithm = "hc",
                         algorithm.args = list(score="bic-g")) #gs, iamb, mmpc and hc

head(str.diff)
attr(str.diff, "threshold")
plot(str.diff)

avg.diff = averaged.network(str.diff)
strength.plot(avg.diff, str.diff, shape = "ellipse") # threshold = attr(str.diff, "threshold")

avg.simpler = averaged.network(str.diff, threshold = 0.25) # change threshold
strength.plot(avg.simpler, str.diff, shape = "ellipse")
