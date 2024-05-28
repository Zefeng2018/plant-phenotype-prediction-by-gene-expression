## load common HVGs
setwd("/home/wuzefeng/MyResearch/Tissue_prediction/")
library(stringr)
orthos <- read.table("orthologs/maize2rice.orthologs.txt",header = T)
orthos$subject_id <- str_replace(string = orthos$subject_id,pattern ="t" ,replacement = "g") 

maize_hvg <- read.table("maize/HVG.txt",header = T)
rice_hvg <- read.table("rice/HVG.txt",header = F)

orthos$maize.hvg <- ifelse(orthos$query_id%in%maize_hvg$Gene,1,0) 
orthos$os.hvg <- ifelse(orthos$subject_id%in%rice_hvg$V1,1,0) 

#common_hvgs <- table(orthos[,c("maize.hvg","os.hvg")])[2,2]
ortho_hvgs <- orthos[orthos$maize.hvg==1&orthos$os.hvg==1,]
ortho_hvgs$orth_hvgs <- paste("ohvg",seq(1:nrow(ortho_hvgs)),sep = "_")


## maize train
maize_HVG_tpm <- read.table("maize/HVG.tpm.txt",header = T)
maize_HVG_tpm <- maize_HVG_tpm[,colnames(maize_HVG_tpm)%in%ortho_hvgs$query_id] # 486 column hvgs
colnames(maize_HVG_tpm) <- ortho_hvgs$orth_hvgs[match(colnames(maize_HVG_tpm),ortho_hvgs$query_id)]  # rename hvgs to ortholog name

maize_sample_info <- read.table("maize/sample.info.txt",sep="\t",header = F)
maize_HVG_tpm$class <- maize_sample_info$V5[match(rownames(maize_HVG_tpm),maize_sample_info$V1)]  # sample * genes * class 


head(sort(table(maize_HVG_tpm$class),decreasing = T),40)
maize_HVG_tpm_sub <- subset(maize_HVG_tpm,maize_HVG_tpm$class%in%c("leaf","root","kernel","seed","endosperm","shoot","stem","anther","embryo"))


# splitting data into 7:3 for training and testing 
library(caret)
set.seed(1)
train_sub = sample(nrow(maize_HVG_tpm_sub),7/10*nrow(maize_HVG_tpm_sub))
train_data = maize_HVG_tpm_sub[train_sub,]
test_data = maize_HVG_tpm_sub[-train_sub,]
train_data$class = as.factor(train_data$class)
test_data$class = as.factor(test_data$class)
train_control <- trainControl(method = "cv",number = 5) 

## make model 1: random forest
library(randomForest)
system.time(rf_model <- randomForest(class~.,data = train_data,importance=T)) # slow 23 min
rf_pred <- predict(rf_model,test_data[-ncol(test_data)]) # 0.9805

## make model 1： naive Bayes 
nb_model <- naiveBayes(class ~ ., data=train_data)
nb_pred <- predict(nb_model, test_data[,-ncol(test_data)], type="class") # ACC:0.8409
saveRDS(nb_model,"models/tissue_nb_models.rds")
loaded_model <-readRDS("/home/wuzefeng/MyResearch/Tissue_prediction/maize/models/tissue_nb_models.rds")



##make model 3： SVM 
library(e1071)
svm_model <- svm(class~., data=train_data,type = 'C',kernel="linear") # 
svm_pred <- predict(object = svm_model,newdata = test_data[-ncol(test_data)])


## make model 4 XGboost
library("xgboost")

train_data = as.matrix(maize_HVG_tpm_sub[train_sub,-ncol(maize_HVG_tpm_sub)])
label_train = as.numeric(as.factor(maize_HVG_tpm_sub$class[train_sub]))-1
dtrain <- xgb.DMatrix(data = train_data, label = label_train)

test_data = as.matrix(maize_HVG_tpm_sub[-train_sub,-ncol(maize_HVG_tpm_sub)])
label_test <- as.numeric(as.factor(maize_HVG_tpm_sub[-train_sub,]$class))-1
dtest <- xgb.DMatrix(data = test_data, label =label_test)

xgb <- xgboost(data = dtrain,max_depth=6, eta=0.5, objective='multi:softmax',num_class=9, nround=25,verbose = 1)


## make model 5 deep learning
library(keras)
library(tensorflow)
HVG_tpm_sub_keras <- maize_HVG_tpm_sub
HVG_tpm_sub_keras$class <- as.numeric(as.factor(HVG_tpm_sub_keras$class))-1
tissues_num_maps <- unique(data.frame(tissues=maize_HVG_tpm_sub$class,num=HVG_tpm_sub_keras$class))

x_train = HVG_tpm_sub_keras[train_sub,-ncol(HVG_tpm_sub_keras)] %>% as.matrix
y_train = HVG_tpm_sub_keras[train_sub,] %>% pull(class) %>% to_categorical(9)
x_test  = HVG_tpm_sub_keras[-train_sub,-ncol(HVG_tpm_sub_keras)] %>% as.matrix
y_test  = HVG_tpm_sub_keras[-train_sub,]  %>% pull(class) %>% to_categorical(9)

dnn_model = keras_model_sequential()
dnn_model %>% 
  layer_dense(units = 1000, activation = 'relu', input_shape = 486) %>% 
  layer_dense(units = 9, activation = 'softmax')
dnn_model %>% summary

dnn_model %>% compile(
  loss      = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics   = c('accuracy')
)
history = dnn_model %>% fit(
  x = x_train, y = y_train,
  epochs           = 200,
  batch_size       = 20,
  validation_split = 0
)
plot(history)
perf = dnn_model %>% evaluate(x_test, y_test) # keras  # acc: 97.6
print(perf)

y_pred  = dnn_model %>% predict(x_test) %>% k_argmax() %>% as.numeric



## 1.  confusion matrix
svm_cm <- confusionMatrix(svm_pred,test_data$class)  # AC:0.9664
xgb_cm <- confusionMatrix(as.factor(predict(xgb, dtest)),as.factor(label_test)) # 0.9799
rf_cm <- confusionMatrix(rf_pred,test_data$class) # 0.9737
nb_cm <- confusionMatrix(nb_pred,test_data$class) # 0.6783

#dnn_cm <- confusionMatrix(as.factor(y_pred),as.factor(HVG_tpm_sub_keras[-train_sub,] %>% pull(class)))
dnn_pre <- as.factor(unlist(lapply(y_pred,function(x)tissues_num_maps$tissues[tissues_num_maps$num==x])))
dnn_true <- as.factor(maize_HVG_tpm_sub[-train_sub,]  %>% pull(class))
dnn_cm <- confusionMatrix(dnn_pre,dnn_true) # 0.9723

#######################################################################################################3
### rice data load
rice_HVG_tpm <- read.table("rice/HVG.tpm.txt",header = T)
rice_HVG_tpm <- rice_HVG_tpm[,colnames(rice_HVG_tpm)%in%ortho_hvgs$subject_id]
colnames(rice_HVG_tpm) <- ortho_hvgs$orth_hvgs[match(colnames(rice_HVG_tpm),ortho_hvgs$subject_id)]
rice_sample_info <- read.table("rice/sample.info.txt",sep="\t",header = F)
rice_HVG_tpm$class <- rice_sample_info$V5[match(rownames(rice_HVG_tpm),rice_sample_info$V1)]

## subset data for ML
head(sort(table(rice_HVG_tpm$class),decreasing = T),40)
rice_HVG_tpm_sub <- subset(rice_HVG_tpm,rice_HVG_tpm$class%in%c("leaf","root","kernel","seed","endosperm","shoot","stem","anther","embryo"))
rice_HVG_tpm_sub <- rice_HVG_tpm_sub[,colnames(maize_HVG_tpm_sub)]

## rf 
maize_rf_pred_rice <- predict(rf_model,rice_HVG_tpm_sub[,-ncol(rice_HVG_tpm_sub)])
maize2rice_rf_cm <- confusionMatrix(maize_rf_pred_rice,as.factor(rice_HVG_tpm_sub$class)) # 0.7629


##  svm
maize_svm_pred_rice <- predict(svm_model,rice_HVG_tpm_sub[,-ncol(rice_HVG_tpm_sub)])
maize2rice_svm_cm <- confusionMatrix(maize_svm_pred_rice,as.factor(rice_HVG_tpm_sub$class)) # 0.5814

## nb
maize_nb_pred_rice <- predict(nb_model,rice_HVG_tpm_sub[,-ncol(rice_HVG_tpm_sub)])
maize2rice_nb_cm <- confusionMatrix(maize_nb_pred_rice,as.factor(rice_HVG_tpm_sub$class)) # 0.503

# xgboost
rice_data = as.matrix(rice_HVG_tpm_sub[,-ncol(rice_HVG_tpm_sub)])
label_rice <- as.numeric(as.factor(rice_HVG_tpm_sub$class))-1
ddata <- xgb.DMatrix(data = rice_data, label =label_rice)
maize2rice_xgb_cm <- confusionMatrix(as.factor(predict(xgb, ddata)),as.factor(label_rice)) # 0.7346
rownames(maize2rice_xgb_cm$table)<-colnames(xgb_cm$table)<-levels(as.factor(rice_HVG_tpm_sub$class))

## dnn cm
maize_dnn_pred_rice  = dnn_model %>% predict(rice_data) %>% k_argmax() %>% as.numeric
dnn_cm <- confusionMatrix(as.factor(maize_dnn_pred_rice),as.factor(label_rice)) # 0.5911

#### plot cm
plot_confusion <- function(cm){
  as.table(cm) %>%
    as_tibble() %>%
    mutate(Prediction=factor(Prediction),
           Truth = factor(Reference,rev(levels(Prediction)))) %>%
    ggplot(aes(Prediction,Truth,fill=log(n+1,2)))+
    geom_tile()+
    geom_text(aes(label=n))+
    scale_fill_gradientn(colors=rev(hcl.colors(9,"Blues")))+
    coord_fixed()+
    theme_minimal()+
    labs(fill = "Number")
}
library(tidyverse)
plot_confusion(xgb_cm)
