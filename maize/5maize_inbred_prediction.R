###############################################################################
## load data
###############################################################################
library(caret)
library(e1071)

setwd("/home/wuzefeng/MyResearch/Tissue_prediction/maize")
# load HVG expression data
HVG_tpm <- read.table("HVG.tpm.txt",header = T)
sample_info <- read.table("sample.info.txt",sep="\t",header = F)
HVG_tpm$class <- sample_info$V5[match(rownames(HVG_tpm),sample_info$V1)]  # sample * genes * class 
######################################
## inbred line prediction
#####################################

HVG_tpm_sub <- subset(HVG_tpm,HVG_tpm$class=="leaf")
HVG_tpm_sub$class <- sample_info$V2[match(rownames(HVG_tpm_sub),sample_info$V1)]
HVG_tpm_sub <- subset(HVG_tpm_sub,HVG_tpm_sub$class%in%c("B73","Wisconsin","A188","B104","W22","Mo17"))
train_control <- trainControl(method = "cv",number = 5,) 
set.seed(1)
train_sub = sample(nrow(HVG_tpm_sub),7/10*nrow(HVG_tpm_sub))
train_data = HVG_tpm_sub[train_sub,]
test_data = HVG_tpm_sub[-train_sub,]
train_data$class = as.factor(train_data$class)
test_data$class = as.factor(test_data$class)

################## 1. SVM
svm_model <- svm(class~., data=train_data,type = 'C',kernel="linear") # 0.9677
svm_pred <- predict(object = svm_model,newdata = test_data[-ncol(test_data)])
svm_result <- confusionMatrix(svm_pred,test_data$stage)


## make model 2: random forest
library(randomForest)
system.time(rf_model <- randomForest(class~.,data = train_data,importance=T)) # slow 23 min
rf_pred <- predict(rf_model,test_data[-ncol(test_data)]) 


## make model 3： naive Bayes 
nb_model <- naiveBayes(class ~ ., data=train_data)
nb_pred <- predict(nb_model, test_data[,-ncol(test_data)], type="class") # ACC:0.8409
#saveRDS(nb_model,"models/tissue_nb_models.rds")
#loaded_model <-readRDS("/home/wuzefeng/MyResearch/Tissue_prediction/maize/models/tissue_nb_models.rds")


## make model 4 XGboost
library("xgboost")

train_data_xg = as.matrix(HVG_tpm_sub[train_sub,-ncol(HVG_tpm_sub)])
label_train = as.numeric(as.factor(HVG_tpm_sub$class[train_sub]))-1

test_data_xg = as.matrix(HVG_tpm_sub[-train_sub,-ncol(HVG_tpm_sub)])
label_test <- as.numeric(as.factor(HVG_tpm_sub[-train_sub,]$class))-1
dtrain <- xgb.DMatrix(data = train_data_xg, label = label_train)
dtest <- xgb.DMatrix(data = test_data_xg, label =label_test)

xgb <- xgboost(data = dtrain,max_depth=6, eta=0.5, objective='multi:softmax',num_class=6, nround=25,verbose = 1)





## make model 6 deep learning
library(keras)
library(tensorflow)
library(tidyverse)
HVG_tpm_sub_keras <- HVG_tpm_sub
HVG_tpm_sub_keras$class <- as.numeric(as.factor(HVG_tpm_sub_keras$class))-1
tissues_num_maps <- unique(data.frame(tissues=HVG_tpm_sub$class,num=HVG_tpm_sub_keras$class))

x_train = HVG_tpm_sub_keras[train_sub,-ncol(HVG_tpm_sub_keras)] %>% as.matrix
y_train = HVG_tpm_sub_keras[train_sub,] %>% pull(class) %>% to_categorical(6)
x_test  = HVG_tpm_sub_keras[-train_sub,-ncol(HVG_tpm_sub_keras)] %>% as.matrix
y_test  = HVG_tpm_sub_keras[-train_sub,]  %>% pull(class) %>% to_categorical(6)

dnn_model = keras_model_sequential()
dnn_model %>% 
  layer_dense(units = 1000, activation = 'relu', input_shape = 2880) %>% 
  layer_dense(units = 6, activation = 'softmax')
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
perf = dnn_model %>% evaluate(x_test, y_test) # keras  # acc: 0.9828
print(perf) # 0.988

y_pred  = dnn_model %>% predict(x_test) %>% k_argmax() %>% as.numeric


######################################################################
# model evaluation
######################################################################
## 1.  confusion matrix
svm_cm <- confusionMatrix(svm_pred,test_data$class)  # AC:0.9677
xgb_cm <- confusionMatrix(as.factor(predict(xgb, dtest)),as.factor(label_test)) # 0.9957
rf_cm <- confusionMatrix(rf_pred,test_data$class) # 0.9784
nb_cm <- confusionMatrix(nb_pred,test_data$class) # 0.8772

dnn_pre <- as.factor(unlist(lapply(y_pred,function(x)tissues_num_maps$tissues[tissues_num_maps$num==x])))
dnn_true <- as.factor(HVG_tpm_sub[-train_sub,]  %>% pull(class))
dnn_cm <- confusionMatrix(dnn_pre,dnn_true) # 0.988
colnames(dnn_cm$table)<-rownames(dnn_cm$table)<-c("A188","B104", "B73","Mo17","W22","Wisconsin")
colnames(xgb_cm$table)<-rownames(xgb_cm$table)<-c("A188","B104", "B73","Mo17","W22","Wisconsin")
## 2. Confuse matrix and visualization
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



## 
df <- data.frame(pre=predict(xgb, dtest),true=label_test)
## 3. ## calculation
library(pROC)
library(tidyverse)

multiclass.roc(as.numeric(test_data$stage),as.numeric(svm_pred), direction = "<")  #  svm: auc 0.9951
multiclass.roc(label_test,predict(xgb, dtest), direction = "<") # 0.9861
multiclass.roc(as.numeric(test_data$stage),as.numeric(rf_pred), direction = "<") # 0.9934
multiclass.roc(as.numeric(dnn_true),as.numeric(dnn_pre), direction = "<")  # 0.9962
multiclass.roc(as.numeric(test_data$stage),as.numeric(nb_pred), direction = "<") # 0.9867


############## pca

library(mixOmics)
info <- HVG_tpm_sub$class
HVG_tpm_sub <- HVG_tpm_sub[,colnames(HVG_tpm_sub)!="class"]
rownames(HVG_tpm_sub) <- NULL

hvg <- list()
hvg$gene <- HVG_tpm_sub
hvg$genotype <- as.factor(info)


MyResult.pca <- pca(hvg$gene)
plotIndiv(MyResult.pca,group = hvg$genotype,title = "Samples",legend = TRUE,size.axis = 20,size.xlabel = 20,size.ylabel = 20)


####### feature importance

model <- xgb.dump(xgb, with_stats= T)
model[1:10] #This statement prints top 10 nodes of the model
# 获得特征的真实名称
names <- colnames(train_data_xg)
# 计算特征重要性矩阵
importance_matrix <- xgb.importance(names, model = xgb)
# 制图
xgb.plot.importance(importance_matrix[1:10,])
# 在最后一步如果失效可能是因为版本问题,你可以尝试:
barplot(importance_matrix[,1])

### gene expression 

subdata <- train_data[,colnames(train_data)%in%c(importance_matrix$Feature[1:5],"class")]
subdata_long <- gather(data = subdata,key = "Gene",value = "Expression",-class)
ggplot(subdata_long,aes(x=class,y=log(Expression+1,2),fill=Gene))+
  geom_boxplot(outliers = TRUE,
               outlier.colour = "gray",
               position=position_dodge(0.7),width=0.5)+
  scale_fill_npg()+
  theme(legend.position = c(0.9,0.85))+
  ylab("Expression / log2(TPM)")+
  xlab("Maize accessions")
