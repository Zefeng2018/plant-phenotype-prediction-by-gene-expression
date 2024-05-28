###############################################################################
## load data
###############################################################################
library(caret)
library(e1071)

setwd("/home/wuzefeng/MyResearch/Tissue_prediction/rice/")
# load HVG expression data
HVG_tpm <- read.table("HVG.tpm.txt",header = T)
sample_info <- read.table("sample.info.txt",sep="\t",header = F)
HVG_tpm$tissue <- sample_info$V5[match(rownames(HVG_tpm),sample_info$V1)]  # sample * genes * class 
HVG_tpm$stage <- sample_info$V4[match(rownames(HVG_tpm),sample_info$V1)] # 21611
######################################
## stage prediction
#####################################
HVG_tpm_sub <- subset(HVG_tpm,HVG_tpm$tissue=="leaf") # 3387
HVG_tpm_sub <- subset(HVG_tpm,HVG_tpm$stage%in%c("seedling","tillering","booting","heading"))
HVG_tpm_sub <- HVG_tpm_sub[,colnames(HVG_tpm_sub)!="tissue"] # 1095


set.seed(1)
train_sub = sample(nrow(HVG_tpm_sub),7/10*nrow(HVG_tpm_sub))
train_data = HVG_tpm_sub[train_sub,]
test_data = HVG_tpm_sub[-train_sub,]
train_data$stage = as.factor(train_data$stage)
test_data$stage = as.factor(test_data$stage)

################## 1. SVM
svm_model <- svm(stage~., data=train_data,type = 'C',kernel="linear") # radial:0.976; sigmoid:0.987; polynomial:0.889; linear:0.998
svm_pred <- predict(object = svm_model,newdata = test_data[-ncol(test_data)])
svm_result <- confusionMatrix(svm_pred,test_data$stage)



## make model 2: random forest
library(randomForest)
system.time(rf_model <- randomForest(stage~.,data = train_data,importance=T)) # slow 23 min
rf_pred <- predict(rf_model,test_data[-ncol(test_data)]) 


## make model 3： naive Bayes 
nb_model <- naiveBayes(stage ~ ., data=train_data)
nb_pred <- predict(nb_model, test_data[,-ncol(test_data)], type="class") # ACC:0.8409
#saveRDS(nb_model,"models/tissue_nb_models.rds")
#loaded_model <-readRDS("/home/wuzefeng/MyResearch/Tissue_prediction/maize/models/tissue_nb_models.rds")


## make model 4 XGboost
library("xgboost")

HVG_tpm_sub_xg <- HVG_tpm_sub
HVG_tpm_sub_xg$stage <- as.numeric(as.factor(HVG_tpm_sub_xg$stage))-1
tissues_num_maps_xg <- unique(data.frame(tissues=HVG_tpm_sub$stage,num=HVG_tpm_sub_xg$stage))

train_data_xg = as.matrix(HVG_tpm_sub[train_sub,-ncol(HVG_tpm_sub)])
label_train = as.numeric(as.factor(HVG_tpm_sub$stage[train_sub]))-1

test_data_xg = as.matrix(HVG_tpm_sub[-train_sub,-ncol(HVG_tpm_sub)])
label_test <- as.numeric(as.factor(HVG_tpm_sub[-train_sub,]$stage))-1
dtrain <- xgb.DMatrix(data = train_data_xg, label = label_train)
dtest <- xgb.DMatrix(data = test_data_xg, label =label_test)

xgb <- xgboost(data = dtrain,max_depth=6, eta=0.5, objective='multi:softmax',num_class=4, nround=25,verbose = 1)
xgb_pre <- predict(xgb, dtest)

## make model 6 deep learning
library(keras)
library(tensorflow)
library(tidyverse)
HVG_tpm_sub_keras <- HVG_tpm_sub
HVG_tpm_sub_keras$stage <- as.numeric(as.factor(HVG_tpm_sub_keras$stage))-1


x_train = HVG_tpm_sub_keras[train_sub,-ncol(HVG_tpm_sub_keras)] %>% as.matrix
y_train = HVG_tpm_sub_keras[train_sub,] %>% pull(stage) %>% to_categorical(4)
x_test  = HVG_tpm_sub_keras[-train_sub,-ncol(HVG_tpm_sub_keras)] %>% as.matrix
y_test  = HVG_tpm_sub_keras[-train_sub,]  %>% pull(stage) %>% to_categorical(4)

dnn_model = keras_model_sequential()
dnn_model %>% 
  layer_dense(units = 1000, activation = 'relu', input_shape = 3997) %>% 
  layer_dense(units = 4, activation = 'softmax')
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
print(perf) # 0.988

y_pred  = dnn_model %>% predict(x_test) %>% k_argmax() %>% as.numeric


######################################################################
# model evaluation
######################################################################
## 1.  confusion matrix
svm_cm <- confusionMatrix(svm_pred,test_data$stage)  # AC:0.9848
xgb_cm <- confusionMatrix(as.factor(unlist(lapply(xgb_pre,function(x)tissues_num_maps_xg$tissues[tissues_num_maps_xg$num==x]))),as.factor(unlist(lapply(label_test,function(x)tissues_num_maps_xg$tissues[tissues_num_maps_xg$num==x])))) # 0.9635
rf_cm <- confusionMatrix(rf_pred,test_data$stage) # 0.9848
nb_cm <- confusionMatrix(nb_pred,test_data$stage) # 0.7872


dnn_pre <- as.factor(unlist(lapply(y_pred,function(x)tissues_num_maps$tissues[tissues_num_maps$num==x])))
dnn_true <- as.factor(HVG_tpm_sub[-train_sub,]  %>% pull(stage))
dnn_cm <- confusionMatrix(dnn_pre,dnn_true) # 0.988
colnames(dnn_cm$table)<-rownames(dnn_cm$table)<-c("V2-stage","V4-stage","V1-stage","V3-stage")
## 2. Confuse matrix and visulization
plot_confusion <- function(cm){
  as.table(cm) %>%
    as_tibble() %>%
    mutate(Prediction=factor(Prediction),
           Truth = factor(Reference,rev(levels(Prediction)))) %>%
    ggplot(aes(Prediction,Truth,fill=n))+
    geom_tile()+
    geom_text(aes(label=n))+
    scale_fill_gradientn(colors=rev(hcl.colors(4,"Blues")))+
    coord_fixed()+
    theme_minimal()+
    labs(fill = "Number")
}
plot_confusion(svm_cm$table)
plot_confusion(dnn_cm)
plot_confusion(xgb_cm)
plot_confusion(rf_cm)
plot_confusion(nb_cm)
## 3. ## calculation
library(pROC)
library(tidyverse)

multiclass.roc(as.numeric(test_data$stage),as.numeric(svm_pred), direction = "<")  #  svm: auc 0.9951
multiclass.roc(label_test,predict(xgb, dtest), direction = "<") # 0.9861
multiclass.roc(as.numeric(test_data$stage),as.numeric(rf_pred), direction = "<") # 0.9934
multiclass.roc(as.numeric(dnn_true),as.numeric(dnn_pre), direction = "<")  # 0.9962
multiclass.roc(as.numeric(test_data$stage),as.numeric(nb_pred), direction = "<") # 0.9867


## 

library(tidyverse)
import <- importance(rf_model)%>%as.data.frame()%>%arrange(-MeanDecreaseAccuracy)


### gene expression plot
library(ggsci)
subdata <- train_data[,colnames(train_data)%in%c(rownames(import)[1:5],"stage")]
subdata_long <- gather(data = subdata,key = "Gene",value = "Expression",-stage)
subdata_long$stage <- factor(subdata_long$stage,levels = c("seedling", "tillering", "booting", "heading"))
ggplot(subdata_long,aes(x=stage,y=log(Expression+1,2),fill=Gene))+
  geom_boxplot(outliers = TRUE,
               outlier.colour = "gray",
               position=position_dodge(0.7),width=0.5)+
  scale_fill_npg()+
  #theme(legend.position = c(0.9,0.85))+
  ylab("Expression / log2(TPM)")+
  xlab("Development stages")+
  theme(legend.position = c(0.85,0.85))










