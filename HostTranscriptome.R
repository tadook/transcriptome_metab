
# Library
library(DESeq2) # BiocManager::install("DESeq2")
library(ggbio) # BiocManager::install("ggbio")
library(tximport) # BiocManager::install("tximport")
library(pheatmap)
library(stringr)
library(tidyverse)
library(ggplot2)
library(readr)
library(org.Hs.eg.db)
library(data.table)
library(haven)
library(rsample)
library(glmnet)
library(ROCR)
library(randomForest)
library(quantable)

# Import
txi_398 <- readRDS("../rawdata/txi_salmon_n398.rds")
transcript_398_df_mod <- read_csv("../rawdata/metadata_DESeq_n398.csv")
metadata <- read_dta("../rawdata/m35_metadata_n1016_2018march.dta")

# Merge files
ddsTxi_virus_398 <- DESeqDataSetFromTximport(txi_398,
                                             colData = transcript_398_df_mod,
                                             design = ~ 1)

keep_virus_398 <- rowSums(counts(ddsTxi_virus_398)) >= 10 # remove mRNA counts < 10
ddsTxi_virus_398_2 <- ddsTxi_virus_398[keep_virus_398,]
ddsTxi_virus_398_2 # 21136 -> 20486 (20596?)
ddsTxi_virus_398_3 <- DESeq(ddsTxi_virus_398_2)

# Make dataset for selecting significant transcriptome
pre_normalized_counts <- t(counts(ddsTxi_virus_398_3, normalized=TRUE)) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "study_id") %>%
  merge(metadata[,c("study_id","CPAPintubate","inpatient_hfo","IntensiveTreatment","intake_sex","Age_mo","raceethn")], by ="study_id") %>%
  mutate(severity = ifelse(CPAPintubate == 1, 1, ifelse(inpatient_hfo == 1, 1, 0))) %>% 
  mutate(intake_sex = ifelse(intake_sex == 2, 1, 0))

tra_variable <- NULL
for (i in 2:20597){ 
  severe.fit <- glm(severity ~ ., 
                    family = binomial(),
                    pre_normalized_counts[,c(i,20601:20604)]) # covariates:intake_sex, Age_mo, raceethn
  sum <- summary(severe.fit)
  sum_tran <- sum$coefficients[2,] %>% as.data.frame() %>% t()
  rownames(sum_tran) <- colnames(pre_normalized_counts)[i]
  tra_variable <- rbind(tra_variable,sum_tran)
}

tra_var_select <- filter(as.data.frame(tra_variable), tra_variable[,4]<0.05)

# Selecting variable and normalize the dataset

normalized_counts <- t(counts(ddsTxi_virus_398_3, normalized=TRUE)) %>% 
  log1p() %>% # normalization: log(x+1)
  as.data.frame() %>%
  dplyr::select(rownames(tra_var_select)) 
  
# MinMax scaling for Deepinsight
MinMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

for (i in 1:849){
  normalized_counts[,i] <- MinMax(normalized_counts[,i])
}

# Final transcriptome dataset 
normalized_counts <- normalized_counts %>% 
  rownames_to_column(var = "study_id") %>%
  merge(metadata[,c("study_id","CPAPintubate","inpatient_hfo","IntensiveTreatment","intake_sex","Age_mo")], by ="study_id") %>%
  mutate(severity = ifelse(CPAPintubate == 1, 1, ifelse(inpatient_hfo == 1, 1, 0))) %>% 
  mutate(intake_sex = ifelse(intake_sex == 2, 1, 0))

fwrite(normalized_counts,"../transcriptome.csv")

# table <- vst(ddsTxi_virus_398_3)
# 
# trans_data <- t(table@assays@data@listData[[1]]) %>%  # http://adv-r.had.co.nz/S4.html
#   as.data.frame() %>%
#   rownames_to_column(var = "study_id") %>%
#   merge(metadata[,c("study_id","CPAPintubate","inpatient_hfo","IntensiveTreatment","intake_sex","Age_mo")], by ="study_id") %>%
#   mutate(severity = ifelse(CPAPintubate == 1, 1, ifelse(inpatient_hfo == 1, 1, 0))) %>% 
#   mutate(intake_sex = ifelse(intake_sex == 2, 1, 0))
# 
# fwrite(trans_data,"trans_data.csv")


# Lasso regression -- simple prediction
# https://bookdown.org/tpinto_home/Regularisation/lasso-regression.html
## Create train data and test data

# normalized_counts <- fread("../normalized_counts.csv")
trans_data <- normalized_counts 
trans_data_pre <- trans_data %>% dplyr::select(-c("study_id","CPAPintubate","inpatient_hfo","IntensiveTreatment","intake_sex","Age_mo")) # "severity"
set.seed(1)
df_split <- initial_split(trans_data_pre, prop = 0.8) 
lasso_train <- training(df_split) 
lasso_test <- testing(df_split) 

lasso_train$severity <- factor(lasso_train$severity) %>% factor(labels = c("No", "Yes")) # severity
lasso_train_x <- lasso_train[,2:161] %>% data.matrix() # 2:20597 2:858
lasso_train_y <- lasso_train[,"severity"] 

lasso_test$severity <- factor(lasso_test$severity) %>% factor(labels = c("No", "Yes")) # severity
lasso_test_x <- lasso_test[,2:161] %>% data.matrix() # 2:20597 2:858
lasso_test_y <- lasso_test[,"severity"] 

# Create model
## Lasso
set.seed(1)
lasso_fit <- cv.glmnet(lasso_train_x,
                       lasso_train_y,
                       family='binomial',
                       type.measure = "class", # auc  # confirm AUC on CV by summary(lasso_fit$cvm)
                       nfolds = 10)
lasso_fit
plot(lasso_fit)
plot(lasso_fit$glmnet.fit, "lambda", label=F)
lasso_fit$lambda.min

lasso_fit <- glmnet(lasso_train_x,
                    lasso_train_y,
                    family='binomial',
                    type.measure = "auc", # auc  # confirm AUC on CV by summary(lasso_fit$cvm)
                    lambda = 1e-12) 

# Check used variables
# df_mapping <- fread("rawdata/metab_dict.csv")
# lasso_val <- coef(lasso_fit, s = "lambda.min")[-1,] %>% as.data.frame()
# colnames(lasso_val) <- "coef"
# lasso_val$metab_id <- rownames(lasso_val)
# lasso_val <- lasso_val %>% 
#   inner_join(df_mapping,by="metab_id") %>% 
#   arrange(desc(abs(coef))) %>% select(metab_id,metab_name,coef)
# fwrite(lasso_val, "lasso_val.csv")

# calculating ROC
score_logit_prob <- predict(lasso_fit,probability=T,lasso_test_x)
score_logit_pred <- prediction(score_logit_prob,lasso_test_y)
rocObj_score_logit <- performance(score_logit_pred,measure="tpr", x.measure="fpr")
prObj_score_logit <- performance(score_logit_pred,measure="prec", x.measure="rec")
aucObj_score_logit <- performance(score_logit_pred,measure="auc")
praucObj_score_logit <- performance(score_logit_pred,measure="aucpr")
auc_score_logit <- aucObj_score_logit@y.values[[1]]
prauc_score_logit <- praucObj_score_logit@y.values[[1]]
auc_score_logit
prauc_score_logit

# ROC drawing
plot(rocObj_score_logit,main="Sick fat prediction: ROC curve", col="red")
# PR drawing
plot(prObj_score_logit,main="Sick fat prediction: PR curve", col="blue",ylim=c(0,1))


## Random Forest
set.seed(1)
rf_fit <- randomForest(severity ~., data=lasso_train, proximity=T) # Error: protect(): protection stack overflow

# calculating ROC
score_logit_prob <- predict(rf_fit,type="prob",lasso_test)
score_logit_pred <- prediction(score_logit_prob[,2],lasso_test$severity)
rocObj_score_logit <- performance(score_logit_pred,measure="tpr", x.measure="fpr")
prObj_score_logit <- performance(score_logit_pred,measure="prec", x.measure="rec")
aucObj_score_logit <- performance(score_logit_pred,measure="auc")
praucObj_score_logit <- performance(score_logit_pred,measure="aucpr")
auc_score_logit <- aucObj_score_logit@y.values[[1]]
prauc_score_logit <- praucObj_score_logit@y.values[[1]]
auc_score_logit
prauc_score_logit

# ROC drawing
plot(rocObj_score_logit,main="Sick fat prediction: ROC curve", col="red")
# PR drawing
plot(prObj_score_logit,main="Sick fat prediction: PR curve", col="blue",ylim=c(0,1))


var_imp <- importance(rf_fit) %>% 
           as.data.frame() %>% 
           rownames_to_column(var = "transcriptome") 
var <-  variable %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "transcriptome") 

var_imp <- merge(var_imp,var,by="transcriptome")
