
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

# df_metabolome_norm <- df_metabolome_combat # for integrated analysis

# extract uid names -------------------------------------------------------

uid <- rownames(df_metabolome_norm)
df_metabolome_norm_uid <- data.frame(uid,df_metabolome_norm)
rownames(df_metabolome_norm_uid) <- NULL

# Import metadata and outcome (asthma) data -------------------------------

df_metadata <- readstata13::read.dta13("../rawdata/m35_metadata_modified2019.dta") %>%
  select(uid, study_id, 
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

transcriptome_variables <- read_csv("../normalized_counts.csv")
metab_data <- metabo_variables %>%  dplyr::select("study_id",rownames(var_select)) 
met_tra_variables <- merge(metab_data, transcriptome_variables, by = "study_id")

fwrite(met_tra_variables,"../normalized_counts.csv")

