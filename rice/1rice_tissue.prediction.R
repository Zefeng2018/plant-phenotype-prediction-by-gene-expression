
library(caret)
library(e1071)

# load data
HVG_tpm <- read.table("HVG.tpm.txt",header = T)
sample_info <- read.table("sample.info.txt",sep="\t",header = F)
HVG_tpm$class <- sample_info$V5[match(rownames(HVG_tpm),sample_info$V1)]

## subset data for ML
head(sort(table(HVG_tpm$class),decreasing = T),40)
HVG_tpm_sub <- subset(HVG_tpm,HVG_tpm$class%in%c("leaf","root","kernel","seed","endosperm","shoot","stem","anther","embryo","spikelet","seedling","panicle"))

## 5 fold cvs
train_control <- trainControl(method = "cv",number = 5,) 

##将数据集分为训练集和测试集,比例为7:3
set.seed(1)
train_sub = sample(nrow(HVG_tpm_sub),7/10*nrow(HVG_tpm_sub))
train_data = HVG_tpm_sub[train_sub,]
test_data = HVG_tpm_sub[-train_sub,]

train_data$class = as.factor(train_data$class)
test_data$class = as.factor(test_data$class)

## make model 2: random forest
library(randomForest)
system.time(rf_model <- randomForest(class~.,data = train_data,importance=T)) # slow 23 min
saveRDS(rf_model,"models/tissue_rf_models.rds")
rf_pred <- predict(rf_model,test_data[-ncol(test_data)]) # 0.9805


## make model 3： naive Bayes 
nb_model <- naiveBayes(class ~ ., data=train_data)
nb_pred <- predict(nb_model, test_data[,-ncol(test_data)], type="class") # ACC:0.8409
saveRDS(nb_model,"tissue_nb_models.rds")
loaded_model <-readRDS("/home/wuzefeng/MyResearch/Tissue_prediction/rice/models/tissue_nb_models.rds")



##make model 4： SVM 
library(e1071)
svm_model <- svm(class~., data=train_data,type = 'C',kernel="linear") # radial:0.976; sigmoid:0.987; polynomial:0.889; linear:0.998
svm_pred <- predict(object = svm_model,newdata = test_data[-ncol(test_data)])
svm_result <- confusionMatrix(svm_pred,test_data$class)
saveRDS(svm_model,"models/tissue_svm_models.rds")


svm_model <- svm(class~., data=train_data,type = 'C',kernel="linear",probability=T) # radial:0.976; sigmoid:0.987; polynomial:0.889; linear:0.998
svm_pred_prob <- predict(object = svm_model,newdata = test_data[-ncol(test_data)],probability = TRUE)
svm_pred_prob <-attr(svm_pred_prob, "probabilities")[,1]
sscurves <- evalmod(scores = svm_pred_prob, labels = ifelse(test_data$class=="leaf",1,-1))
autoplot(sscurves)

unknown_data <- subset(HVG_tpm,HVG_tpm$class=="/")
svm_unknown_pred <- predict(object = svm_model,newdata = unknown_data[-ncol(unknown_data)])


## make model 5 XGboost
library("xgboost")

HVG_tpm_sub_xg <- HVG_tpm_sub
HVG_tpm_sub_xg$class <- as.numeric(as.factor(HVG_tpm_sub_xg$class))-1
tissues_num_maps_xg <- unique(data.frame(tissues=HVG_tpm_sub$class,num=HVG_tpm_sub_xg$class))


train_data = as.matrix(HVG_tpm_sub[train_sub,-ncol(HVG_tpm_sub)])
label_train = as.numeric(as.factor(HVG_tpm_sub$class[train_sub]))-1

test_data = as.matrix(HVG_tpm_sub[-train_sub,-ncol(HVG_tpm_sub)])
label_test <- as.numeric(as.factor(HVG_tpm_sub[-train_sub,]$class))-1

dtrain <- xgb.DMatrix(data = train_data, label = label_train)
dtest <- xgb.DMatrix(data = test_data, label =label_test)

xgb <- xgboost(data = dtrain,max_depth=6, eta=0.5, objective='multi:softmax',num_class=12, nround=25,verbose = 1)
xgb_pre <- predict(xgb, dtest)

## make model 6 deep learning
library(keras)
library(tensorflow)
HVG_tpm_sub_keras <- HVG_tpm_sub
HVG_tpm_sub_keras$class <- as.numeric(as.factor(HVG_tpm_sub_keras$class))-1
tissues_num_maps <- unique(data.frame(tissues=HVG_tpm_sub$class,num=HVG_tpm_sub_keras$class))

x_train = HVG_tpm_sub_keras[train_sub,-ncol(HVG_tpm_sub_keras)] %>% as.matrix
y_train = HVG_tpm_sub_keras[train_sub,] %>% pull(class) %>% to_categorical(12)
x_test  = HVG_tpm_sub_keras[-train_sub,-ncol(HVG_tpm_sub_keras)] %>% as.matrix
y_test  = HVG_tpm_sub_keras[-train_sub,]  %>% pull(class) %>% to_categorical(12)

dnn_model = keras_model_sequential()
dnn_model %>% 
  layer_dense(units = 1000, activation = 'relu', input_shape = 3997) %>% 
  layer_dense(units = 12, activation = 'softmax')
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
#################################################
# feature importance
#################################################
plot(varImp(rpart_model))
varImpPlot(rf_model)

######################################################################
# model evaluation
######################################################################
## 1.  confusion matrix
     
svm_cm <- confusionMatrix(svm_pred,test_data$class)  # AC:0.9585
xgb_cm <- confusionMatrix(as.factor(unlist(lapply(xgb_pre,function(x)tissues_num_maps_xg$tissues[tissues_num_maps_xg$num==x]))),as.factor(unlist(lapply(label_test,function(x)tissues_num_maps_xg$tissues[tissues_num_maps_xg$num==x])))) # 0.9675
rf_cm <- confusionMatrix(rf_pred,test_data$class) # 0.9551
nb_cm <- confusionMatrix(nb_pred,test_data$class) # 0.5939

#dnn_cm <- confusionMatrix(as.factor(y_pred),as.factor(HVG_tpm_sub_keras[-train_sub,] %>% pull(class)))
dnn_pre <- as.factor(unlist(lapply(y_pred,function(x)tissues_num_maps$tissues[tissues_num_maps$num==x])))
dnn_true <- as.factor(HVG_tpm_sub[-train_sub,]  %>% pull(class))
dnn_cm <- confusionMatrix(dnn_pre,dnn_true) # 0.0.9608

## 2. Confuse matrix and visulization
plot_confusion <- function(cm){
  as.table(cm) %>%
    as_tibble() %>%
    mutate(Prediction=factor(Prediction),
           Truth = factor(Reference,rev(levels(Prediction)))) %>%
    ggplot(aes(Prediction,Truth,fill=log(n+1,2)))+
    geom_tile()+
    geom_text(aes(label=n))+
    scale_fill_gradientn(colors=rev(hcl.colors(12,"Blues")),
                         breaks=seq(0,log(1018+1,2),2),
                         labels=c(0,250,500,750,1000))+
    coord_fixed()+
    theme_minimal()+
    labs(fill = "Number")
}

plot_confusion(rpart_cm)
plot_confusion(svm_cm$table)
plot_confusion(dnn_cm)
plot_confusion(xgb_cm)
plot_confusion(rf_cm)
plot_confusion(nb_cm)

## 3. ## calculation
library(pROC)
library(tidyverse)

multiclass.roc(as.numeric(test_data$class),as.numeric(svm_pred), direction = "<")  #  svm: auc 0.9717
multiclass.roc(HVG_tpm_sub_keras[-train_sub,] %>% pull(class),y_pred, direction = "<") # dnn: 0.9602
multiclass.roc(label_test,predict(xgb, dtest), direction = "<") # 0.9777
multiclass.roc(as.numeric(test_data$class),as.numeric(rf_pred), direction = "<") # 0.9442
multiclass.roc(as.numeric(test_data$class),as.numeric(nb_pred), direction = "<") # 0.8914







