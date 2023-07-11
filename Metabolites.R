
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(dplyr)
library(robustbase)
library(StatMatch)
library(sva) # BiocManager::install("sva")
library(RColorBrewer)
library(colourvalues)
library(readxl)
library(readr)
library(Rcpm) # devtools::install_github("ricoderks/Rcpm") #Git needs to be installed
library(bnlearn)
library(Rgraphviz) # BiocManager::install("Rgraphviz")
library(splitTools) 
library(ranger)
library(randomForest)
library(ROCR)
library(data.table)
library(rsample)

# Import dataset ----------------------------------------------------------

df_metabolome <- read_csv("../rawdata/metabodata_npa_metab_id.csv") %>% arrange(uid) 


# Batch effect corection --------------------------------------------------

df_mapping <- read_csv("../rawdata/mapping_metab_id_name.csv") 

covadjdata <- read_xlsx("../rawdata/m35_metadata_n1016_2023.03.28.xlsx")  
batchdata <- read.csv(file="../rawdata/batch_npa.csv", header=TRUE, sep=",",check.names=FALSE)
phenodata <-  batchdata %>% merge(covadjdata,by="study_id",all.y=F)
phenodata1 <- phenodata[,-2]
rownames(phenodata1) <- phenodata[,2]

normdata <- read.csv(file="../rawdata/metabodata_npa_metab_id_forcombat_n1013.csv", row.names=1, header=TRUE)
normdatamat <- as.matrix(normdata)
batch = phenodata1$batch
modcombat = model.matrix(~IntensiveTreatment+RVonly_NPA, data=phenodata1) # "mod" include 1 or 2 variables associated with batch number.
combat_normdata1 = ComBat(dat=normdatamat, mod=modcombat, batch=batch, par.prior=TRUE, prior.plots=FALSE) # mod=modcombat
df_metabolome_combat<-t(combat_normdata1)  # data includes negative values. This limits available normalization methods.


# Normalization for MWAS------------------------------------------------------

## PQN normalization
df_metabolome_norm <- pqn(df_metabolome_combat)

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

# Normalization for Deepinsight--------------------------------------------

## Merge study_id
data_metabolome <- metabo_variables %>%  dplyr::select("study_id",rownames(var_select)) 
data_metabolome_norm <- data_metabolome[-1]

## Check the distribution after normalization
hist(data_metabolome_norm[,4])

## Robust scaling for the Deepinsight analysis
Robust_scale <- function(x) {
  interquartile = quantile(x)[4] - quantile(x)[2]
  (x - median(x)) / interquartile
}

for (i in 1:128){
  data_metabolome_norm[,i] <- Robust_scale(data_metabolome_norm[,i])
}

## Convert the outlier
Conv_outlier <- function(x) {
  interquartile = quantile(x)[4] - quantile(x)[2]
  quartile_min = quantile(x)[2] - 1.5 * interquartile
  quartile_max = quantile(x)[4] + 1.5 * interquartile
  x[ x < quartile_min ] <- quartile_min
  x[ x > quartile_max ] <- quartile_max
  x
}

for (i in 1:128){
  data_metabolome_norm[,i] <- Conv_outlier(data_metabolome_norm[,i])
}

## Min-max scale
MinMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

for (i in 1:128){
  data_metabolome_norm[,i] <- MinMax(data_metabolome_norm[,i])
}

## Check the distribution after normalization
hist(data_metabolome_norm[,4])

## Merge with study_id
metab_data <- cbind(data_metabolome[1],data_metabolome_norm)

# Check the missing value
table(is.na(metab_data))

for (i in 2:129){
  names(metab_data)[which(names(metab_data)==colnames(metab_data[i]))] <- df_mapping[charmatch(colnames(metab_data[i]),df_mapping$id_seq_metab),2]
}

fwrite(metab_data,"../metabolites.csv")

# Merge metabolites and transcriptome data --------------------------------

metab_data <- read_csv("../metabolites.csv")
transcriptome_variables <- read_csv("../transcriptome.csv")

metab_trans <- merge(metab_data, transcriptome_variables, by = "study_id")
metab_trans <- metab_trans[,c(1,718:723,2:717)]
pre_metab_trans <- metab_trans

table(is.na(metab_trans))

fwrite(metab_trans,"../trans_metab.csv")

# Make interaction term ---------------------------------------------------

for (i in 8:135){ # This loop takes about 15 minutes
  for (k in 136:723){
    add_column <- data.frame(metab_trans[,i] * metab_trans[,k])
    colnames(add_column) <- paste0(colnames(metab_trans)[i],"-",colnames(metab_trans)[k])
    metab_trans <- cbind(metab_trans, add_column)
  }
} 

## Selecting significant interaction terms

metab_trans_sel <- merge(df_metadata[,c("study_id","raceethn")], metab_trans, by = "study_id")

variable_mettra <- NULL
for (i in 9:75988){ 
  severe.fit <- glm(severity ~ ., 
                    family = binomial(),
                    metab_trans_sel[,c(i,2,6:8)]) # covariates:intake_sex, Age_mo, raceethn
  sum <- summary(severe.fit)
  sum_mettra <- sum$coefficients[2,] %>% as.data.frame() %>% t()
  rownames(sum_mettra) <- colnames(metab_trans_sel)[i]
  variable_mettra <- rbind(variable_mettra,sum_mettra)
}

variable_mettra <- as.data.frame(variable_mettra)[717:75980,]
fdr <- p.adjust(variable_mettra[,4], method ="BH")
var_mettra_select <- cbind(as.data.frame(variable_mettra), as.data.frame(fdr))
var_mettra_select <- filter(var_mettra_select, var_mettra_select[,5]<0.01)

metab_trans <- metab_trans %>% dplyr::select(colnames(metab_trans[1:723]),rownames(var_mettra_select)) 

## Normalization for interaction terms

## Check the distribution before normalization
hist(metab_trans[,1000])

## Robust scaling for the Deepinsight analysis
Robust_scale <- function(x) {
  interquartile = quantile(x)[4] - quantile(x)[2]
  x[ quantile(x)[4] == 0 ] <- x  # for zero skewed variable
  x[ quantile(x)[4] != 0 ] <- (x - median(x)) / interquartile
  x
}

for (i in 724:5558){
  metab_trans[,i] <- Robust_scale(metab_trans[,i])
}

## Convert the outliers
Conv_outlier <- function(x) {
  interquartile = quantile(x)[4] - quantile(x)[2]
  quartile_min = quantile(x)[2] - 1.5 * interquartile
  quartile_max = quantile(x)[4] + 1.5 * interquartile
  x[ quantile(x)[4] != 0 & x < quartile_min ] <- quartile_min
  x[ quantile(x)[4] != 0 & x > quartile_max ] <- quartile_max
  x
}

for (i in 724:5558){
  metab_trans[,i] <- Conv_outlier(metab_trans[,i])
}

## Min-max scale
MinMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

for (i in 724:5558){
  metab_trans[,i] <- MinMax(metab_trans[,i])
}

## Check the distribution after normalization
hist(metab_trans[,1000])

# Check the missing value
table(is.na(metab_trans))

## Save the data 
metab_trans_int <- metab_trans
fwrite(metab_trans_int,"../trans_metab_int.csv")


# Check the precision of RF -----------------------------------------------

metab_trans <-  fread("../trans_metab_int.csv")
  
# metab_trans <- merge(metab_data, transcriptome_variables, by = "study_id")
# metab_trans <- metab_trans[,c(1,979:984,2:978)]
# metab_trans <- metab_trans[,-1:-6]
# vcol <- paste(rep("v",977), seq(1,977,1), sep = "", collapse = ",")
# vcol <- str_split(vcol,",")
# vcol <- vcol[[1]]
# colnames(metab_trans) <- c("severity", vcol)
# metab_trans$severity <- as.factor(metab_trans$severity)

metab_trans <- metab_trans[,-1:-6]
vcol <- paste(rep("v",5551), seq(1,5551,1), sep = "", collapse = ",")
vcol <- str_split(vcol,",")
vcol <- vcol[[1]]
colnames(metab_trans) <- c("severity", vcol)
metab_trans$severity <- as.factor(metab_trans$severity)

# Split data into partitions
set.seed(2)
split_strat <- initial_split(metab_trans, prop = 0.8,
                             strata = 'severity')
train_data <- training(split_strat)
test_metab_trans <- testing(split_strat)

split_strat2 <- initial_split(test_metab_trans, prop = 0.5,
                             strata = 'severity')
val_data <- training(split_strat2)
test_data <- testing(split_strat2)

table(train_data$severity)
table(val_data$severity)
table(test_data$severity)

## Import dataset for metab_trans_int
X_train <- fread("outputs/X_train6_data.csv")
X_val <- fread("outputs/X_val2_data.csv")
X_test <- fread("outputs/X_test2_data.csv")
y_train <- fread("outputs/y_train6_data.csv")
y_val <- fread("outputs/y_val2_data.csv")
y_test <- fread("outputs/y_test2_data.csv")

train_data <- cbind(y_train[,-1], X_train[,-1])
val_data <- cbind(y_val[,-1], X_val[,-1])
test_data <- cbind(y_test[,-1], X_test[,-1])

colnames(train_data) <- c("severity", vcol)
colnames(val_data) <- c("severity", vcol)
colnames(test_data) <- c("severity", vcol)
train_data$severity <- as.factor(train_data$severity)
val_data$severity <- as.factor(val_data$severity)
test_data$severity <- as.factor(test_data$severity)

## Random Forest
set.seed(1)
rf_fit <- randomForest(severity ~., data=train_data, proximity=T) 

## calculating AUC
### validation data
score_prob <- predict(rf_fit,type="prob",val_data)
score_pred <- prediction(score_prob[,2],val_data$severity)
rocObj_score <- performance(score_pred,measure="tpr", x.measure="fpr")
aucObj_score <- performance(score_pred,measure="auc")
auc_score <- aucObj_score@y.values[[1]]
auc_score

### test data
score_prob <- predict(rf_fit,type="prob",test_data)
score_pred <- prediction(score_prob[,2],test_data$severity)
rocObj_score <- performance(score_pred,measure="tpr", x.measure="fpr")
aucObj_score <- performance(score_pred,measure="auc")
auc_score <- aucObj_score@y.values[[1]]
auc_score



# Sensitivity analyses ---------------------------------------------

metab_trans <-  fread("../trans_metab_int.csv")
feature_severe <- read.table("../feature_severe.csv", header = FALSE, sep = ",")

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

severe_variables <- as.vector(feature_severe)
metab_trans_severe <- select(metab_trans, c("study_id",severe_variables[[1]])) %>% 
                      merge(df_metadata[c("study_id","RSV","HRV","Age_mo","intake_sex","raceethn")], by = "study_id") %>% 
                      as.data.frame()

metab_trans_severe$RSV <- as.factor(metab_trans_severe$RSV)
metab_trans_severe$HRV <- as.factor(metab_trans_severe$HRV)
metab_trans_severe$intake_sex <- as.factor(metab_trans_severe$intake_sex)
metab_trans_severe$raceethn <- as.factor(metab_trans_severe$raceethn)

## RSV
RSV_mettra <- NULL
for (i in colnames(metab_trans_severe[2:765])) { 
  formula_str <- paste("`", i, "` ~ RSV + Age_mo + intake_sex + raceethn", sep = "")
  RSV.fit <- lm(as.formula(formula_str), 
                data = metab_trans_severe[,c(i,"RSV","Age_mo","intake_sex","raceethn")])
  sum <- summary(RSV.fit)
  sum_mettra <- sum$coefficients[2,] %>% as.data.frame() %>% t()
  rownames(sum_mettra) <- i
  RSV_mettra <- rbind(RSV_mettra,sum_mettra)
}

fdr <- p.adjust(RSV_mettra[,4], method ="BH")
RSV_mettra_select <- cbind(as.data.frame(RSV_mettra), as.data.frame(fdr))
RSV_mettra_select <- filter(RSV_mettra_select, RSV_mettra_select[,4] < 0.1)

write_csv(as.data.frame(rownames(RSV_mettra_select)),"RSV_mettra_01.csv")

## HRV
HRV_mettra <- NULL
for (i in colnames(metab_trans_severe[2:765])) { 
  formula_str <- paste("`", i, "` ~ HRV + Age_mo + intake_sex + raceethn", sep = "")
  HRV.fit <- lm(as.formula(formula_str), 
                data = metab_trans_severe[,c(i,"HRV","Age_mo","intake_sex","raceethn")])
  sum <- summary(HRV.fit)
  sum_mettra <- sum$coefficients[2,] %>% as.data.frame() %>% t()
  rownames(sum_mettra) <- i
  HRV_mettra <- rbind(HRV_mettra,sum_mettra)
}

fdr <- p.adjust(HRV_mettra[,4], method ="BH")
HRV_mettra_select <- cbind(as.data.frame(HRV_mettra), as.data.frame(fdr))
HRV_mettra_select <- filter(HRV_mettra_select, HRV_mettra_select[,4] < 0.1)

write_csv(as.data.frame(rownames(HRV_mettra_select)),"HRV_mettra_01.csv")
