library(data.table)
library(dplyr)
library(smotefamily)

X_train6_data <- fread("outputs/X_train6_data.csv")
y_train6_data <- fread("outputs/y_train6_data.csv")
X_train8_data <- fread("outputs/X_train8_data.csv")
y_train8_data <- fread("outputs/y_train8_data.csv")

train6_data <- merge(X_train6_data, y_train6_data, by = "V1") %>% select(-V1)
train8_data <- merge(X_train8_data, y_train8_data, by = "V1") %>% select(-V1)

table(train6_data$severity)

smote_train6_data <- SMOTE(train6_data, train6_data$severity, dup_size = 50) 
smote_train6_data <- data.table(smote_train6_data$data) %>% select(-class)
smote_train6_data <- SMOTE(smote_train6_data, smote_train6_data$severity, dup_size = 0) 
table(smote_train6_data$data$severity)
smote_train6_data <- data.table(smote_train6_data$data) %>% select(-class)

smote_train8_data <- SMOTE(train8_data, train8_data$severity, dup_size = 50) 
smote_train8_data <- data.table(smote_train8_data$data) %>% select(-class)
smote_train8_data <- SMOTE(smote_train8_data, smote_train8_data$severity, dup_size = 0) 
table(smote_train8_data$data$severity)
smote_train8_data <- data.table(smote_train8_data$data) %>% select(-class)

X_train6_new <- smote_train6_data[,-"severity"]
y_train6_new <- smote_train6_data[,"severity"]
X_train8_new <- smote_train8_data[,-"severity"]
y_train8_new <- smote_train8_data[,"severity"]

fwrite(X_train6_new,"outputs/X_train6_smote.csv", row.names = T)
fwrite(y_train6_new,"outputs/y_train6_smote.csv", row.names = T)
fwrite(X_train8_new,"outputs/X_train8_smote.csv", row.names = T)
fwrite(y_train8_new,"outputs/y_train8_smote.csv", row.names = T)
