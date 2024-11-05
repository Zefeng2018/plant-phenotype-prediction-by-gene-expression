# Plant phenotype prediction by gene expression

## Inroduction
Research on the dynamic expression of genes in plants is important for understanding different biological processes. In this study, we applied ML approach to a large-scale well-processed gene expression dataset of the two crops, maize and rice, to model the relationship between the gene expression patterns and phenotype (e.g. tissue types, development stages, stress types, etc.). The figure below described the workflow.

![image](https://github.com/Zefeng2018/plant-phenotype-prediction-by-gene-expression/blob/main/images/img.png)

## Script usage

This pipeline was mainly performed by R programming (R version 4.3). There are two folders in this repository, the maize and rice folders. In maize folder, the main functions of each script are as followings.

#### 0HVG_identification.R
This R script is used to identify HGVs by taking maize gene expression matrix (expression levels of genes across different samples) as the input. The output file is a data frame that contains the detailed information about each HGV. 
#### 1maize_tissue_prediction.R
This R script is used to model the relationship between the gene expression level of HGVs and tissue types.
#### 2maize_tissue_prediction_shuffling.R
This script is used to model the relationship between the gene expression level of HVGs and shuffled tissue types.
#### 3maize_tissue_prediction_random_sampling.R
This script is used to model the relationship between the gene expression level of HVGs and shuffled tissue types.
#### 4maize_development_stage.prediction.R
This script is used to model the relationship between the gene expression level of HVGs and maize development stages.
#### 5maize_inbred_prediction.R
This script is used to model the relationship between the gene expression level of HVGs and maize inbred lines or cultivars.
#### 6maize_stress_prediction.R
This script is used to model the relationship between the gene expression level of HVGs and the type of stress experienced by maize
#### 7maize_rice_HVG_compare.R
This script is used to determine the conservation of HVGs between maize and rice by random sampling.
#### 8.maize_rice_model_compare.R
This script is used to predict rice tissue types using the models trained from maize dataset.
#### 9rice_maize_model_compare.R
This script is used to predict maize tissue types using the models trained from rice dataset.

In the rice folder, the main functions of each script are as follows.
#### 0HVG_identification

#### 1rice_tissue.prediction

#### 2rice_development_stage.prediction

#### 3rice_inbred_prediction

#### 4stress_prediction
