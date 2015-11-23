---
layout: projects
authors: Matahi MOARII, Fabien REYAL, Jean-Philippe VERT 
title: Epigenomic alterations in breast carcinoma from primary tumor to locoregional recurrences 
location: Mines Paristech, Institut Curie.
categories: projects-2014
---

# Introduction

This page contains the pipeline analysis for the article [Epigenomic alterations in breast carcinoma from primary tumor to locoregional recurrences](link.zip). Source code is available [here](code/CpG.zip). The whole analysis can be computed by running the **run_all.R**

## Data processing

## Data analysis 

### Part A: Analysis of average CGI+SS patterns

In the first part, we compare the average CGI+SS profiles for a given dataset (**i.e** Cancerous breast) across the genome to assess whether specific CGI+SS profiles exist and whether they are associated with any specific gene expression levels.

```{r}
DiseaseList <- c('BRCA','COAD','LUAD')
Type <- c('Cancerous','Normal')
```


**A.1)** Filter CGI+SS with at least 20 probes:

```{r}
load("../../data/processed/fData/CpGIslands_probe_size.RData")
list_big_island <- which(CpGIslands.probesize >=20)
```
This reduces the number of CGIs studied from 27K to 1827 CGIs.

**A.2)** For each type of tissue and each CGI+SS, we calculate a probewise average profile:

```{r}
source('fun/calculate_Mean_PC.R')
for (DiseaseName in DiseaseList)
{
        out <- calculate_Mean_PC(Disease=DiseaseName,type="Cancerous",proc=F)
        out <- calculate_Mean_PC(Disease=DiseaseName,type="Normal",proc=F)
}
```

**A.3)** We then perform dynamic time warping to assess for a given tissue and type (normal or cancerous) the distance between two different CGI+SS profiles:

```{r}
source('fun/calculate_Mean_PC.R')
for (DiseaseName in DiseaseList)
{
        out <- calculate_Mean_PC(Disease=DiseaseName,type="Cancerous",proc=F)
        out <- calculate_Mean_PC(Disease=DiseaseName,type="Normal",proc=F)
}
```

This outputs a 1827 x 1827 matrix that gives the DTW distance between all CGI+SS profiles

**A.4)** We then perform a hierarchical clustering (linkage=Ward):

```{r}
source('fun/CGI_analysis.R')
for (DiseaseName in DiseaseList)
{
        out <- analyze_CGI_clusters(Disease=DiseaseName,cutoff=3,type="Cancerous")
        out <- analyze_CGI_clusters(Disease=DiseaseName,cutoff=2,type="Normal")
        ## Value for cutoff (i.e number of clusters is given by the hierarchical clustering)
}
```

We observe 2 and 3 clusters of CGI+SS profiles for normal and cancerous tissues respectively.

**A.5)** We plot the characteristic profiles in each cluster i.e the CGI+SS profiles with the lowest mean distance with other CGI+SS profiles in the cluster:

```{r}
source('fun/plot_characteristic_profiles.R')
for (DiseaseName in DiseaseList)
{
        out <- plot_characteristic_profiles(Disease=DiseaseName,type="Cancerous")
        out <- plot_characteristic_profiles(Disease=DiseaseName,type="Normal")
}
```

**A.6)** We assess whether a given CGI+SS is clustered in the same cluster in normal or cancerous tissues:

```{r}
source('fun/compare_clusters.R')
for (DiseaseName in DiseaseList)
{
        compare_clusters(Disease1=DiseaseName,type1="Cancerous", analysis="Mean")
}
```

Clusters are mostly stable between normal and cancerous tissues beside the cancerous-specific cluster that is derived from CGI+SS coming from cluster 1 and 2 in normal tissues.

**A.7)** We then assess whether a given CGI+SS is clustered in the same cluster between different tissues:

```{r}
DiseaseListbis <- c(DiseaseList, DiseaseList[1])
for (k in 1:(length(DiseaseList)-1))
{
        compare_clusters(Disease1=DiseaseList[k],Disease2=DiseaseList[k+1],type1="Normal", analysis="Mean")
        compare_clusters(Disease1=DiseaseList[k],Disease2=DiseaseList[k+1],type1="Cancerous", analysis="Mean")
}
```
Clusters are less stable between tissues!

**A.8)** Finally we look at the link between the CGI+SS patterns and gene expression levels

```{r}
source('fun/compare_GE_clusters.R')
for (DiseaseName in DiseaseList)
{
        compare_GE_clusters(DiseaseName= DiseaseName)
}
```

- In normal tissue, CGI+SS in cluster 2 (hypermethylated CGIs) are not significantly less expressed than CGI+SS in cluster 1 (hypomethylated CGIs)
- In cancerous tissues, we observe that genes associated with cluster 3 are significantly repressed compared to genes associated with cluster 1 and 2.
- Using the clustering of CGI+SS associated with cancerous tissues, and looking at the gene expression distribution in normal tissues, we observe that the genes associated with cluster 3 are also repressed in normal tissues (although the CGI+SS clustering in normal tissues do not have a cluster 3 i.e the CGI+SS are either hypo/hypermethylated). This suggests that overall methylation variations might not be causal in the repression of the genes associated.

### Part B: Inter-individual methylation variations to predict gene expression variations

Average methylation patterns were not associated with gene expression variations. In the second part, we look the power of inter-individual methylation variations (in the CGI+SS) in a specific dataset, to predict the gene expression variations of the associated genes.

**B.1)** We build a regression setting where we assess the predictive power of methylation variations to predict gene expression variations.

We assess the predictive power by performing, for each dataset and for each CGI+SS, a cross-validation procedure (nfolds=3) where we train the parameter of a Lasso on 2/3 of the dataset and we test the prediction on the remaining 1/3 of the dataset.
The performance is assessed with R^2= cor(yhat, ytest)^2 which is a value between 0 and 1 with 1 being the highest. We bootstrap the prediction procedure (nboostrap=100) and we get a final average R^2 for each gene.

We assess the predictive power using only the mean CGI information or the full CGI+SS information. Supplementary analyses include all the CGIs associated with a gene or taking the full methylome to predict the gene expression or just the methylation level of the associated chromosome.

```{r}
source("fun/predict_GE.R")
for (DiseaseName in DiseaseList)
{
        predict_GE(DiseaseName= DiseaseName, type="Cancerous", preprocessing="CGIs", MethylationAnalysis="Mean")
        predict_GE(DiseaseName= DiseaseName, type="Normal", preprocessing="CGIs", MethylationAnalysis="Mean")


        predict_GE(DiseaseName= DiseaseName, type="Cancerous", preprocessing="CGIs", MethylationAnalysis="Promoter")
        predict_GE(DiseaseName= DiseaseName, type="Normal", preprocessing="CGIs", MethylationAnalysis="Promoter")

}
```

**B.2)** Summary of the results:

```{r}
source('fun/analyze_prediction.R')
for (DiseaseName in DiseaseList)
{
        analyze_prediction(DiseaseName)
}
```

**B.3)** We also had the CNV information in the regression model to assess whether the performance in improved (nfolds=3, nboostrap=100):

```{r}
source("fun/predict_GE_CNV.R")
for (DiseaseName in DiseaseList)
{
        predict_GE_CNV(DiseaseName= DiseaseName, type="Cancerous", preprocessing="CGIs", MethylationAnalysis="Mean")
        predict_GE_CNV(DiseaseName= DiseaseName, type="Normal", preprocessing="CGIs", MethylationAnalysis="Mean")


        predict_GE_CNV(DiseaseName= DiseaseName, type="Cancerous", preprocessing="CGIs", MethylationAnalysis="Promoter")
        predict_GE_CNV(DiseaseName= DiseaseName, type="Normal", preprocessing="CGIs", MethylationAnalysis="Promoter")

}
```
**B.4)** Summary of the results:

```{r}
source('fun/analyze_prediction_CNV.R')
for (DiseaseName in DiseaseList)
{
        analyze_prediction_CNV(DiseaseName)
}
```

**B.5)** We then compare the prediction performance with noCNV info:

```{r}
source('fun/compare_prediction_Normal_Cancerous.R')
for (DiseaseName in DiseaseList)
{
        compare_prediction_CNV_noCNV(DiseaseName)
}
```

**B.6)** We compare the prediction performance between different tissues:

```{r}
source('fun/compare_prediction_interCancer.R')
compare_prediction_interCancer()
```

