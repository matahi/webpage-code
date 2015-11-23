---
layout: projects
authors: Matahi MOARII, Fabien REYAL, Jean-Philippe VERT 
title: CIMP analysis 
location: Mines Paristech, Institut Curie.
categories: projects-2015
---

# Introduction

This page contains the pipeline analysis for the "CIMP analysis". Source code is available [here](code/CpG.zip). The whole analysis can be computed by running the **run_all.R** after data processing.

## Data processing

Data processing is done by running the **run_all.R** located in **data/src**. The procedure is standard (see report).

## Data analysis 

### Part A: Analysis of average CGI+SS patterns

```{r}
DiseaseList <- c(('BLCA','BRCA','COAD','LUAD','STAD')
```

**1)** We first assess CIMP in each tissue using the same methodology on genome-wide methylation profiles by performing hierarchical clustering on the top 5% most variant probes in each disease.

```{r}
source('fun/analyze_CIMP_all_CGIs.R')
for (DiseaseName in DiseaseList)
{
        print(DiseaseName)
        out <- analyze_CIMP_all_CGIs(DiseaseName=DiseaseName,CIMP.Number = 2, calc.Var= T)
}
```


**2)** We assess the robustness of the clusters by varying the number of CGIs considered from 1 to 10 percent. At the same time, we also look at the stability of 3 clusters to assess the existence of a CIMP-low phenotype.

```{r}
var.list <- c(1,2,5,10)
CIMP.list <- c(2,3)
source('fun/analyze_CIMP_all_CGIs.R')
for (CIMP.Number in CIMP.list)
{
        for (var.thresh in var.list)
        {
                for (DiseaseName in DiseaseList)
                {
                        print(paste0('Analyzing ',DiseaseName,'...', ' with var=',var.thresh,'% and CIMP.number=',CIMP.Number))
                        out <- analyze_CIMP_all_CGIs(DiseaseName=DiseaseName,CIMP.Number = CIMP.Number, calc.Var= F, var.thresh = var.thresh)
                }
        }
}

### 2.A) Look at cluster robustness
source('fun/cluster_analysis.R')
out <- cluster_analysis(DiseaseList=DiseaseList, var.thresh=var.thresh)

### 2.B) Look at cluster robustness given the var.thresh
source('fun/cluster_analysis_var.R')
var.list <- c(1,2,5,10)

for (Disease in DiseaseList)
{
        out <- cluster_analysis_var(DiseaseName=Disease, CIMP.Number=3, var.list= var.list)
}
```

**3)** We fix the top 5% CGIs to define the CIMP-signature instead of another cutoff as a tradeoff between relevant probes and having a wide enough coverage. We then analyze whether there is a common panel of probes between the tissue-specific CIMP-signature.

```{r}
var.thresh <- 5
source('fun/compare_panel_all_CGIs.R')
DiseaseList <- c('BRCA','BLCA','COAD','LUAD','STAD')
out <- compare_panel_all_CGIs(DiseaseList, var.thresh=var.thresh)
```

We obtain a subset of 89 CGIs common between all the CIMP-signatures.

**4)** By combining the samples from the different tissues, we then perform clustering on this common CIMP-signature:

```{r}
source('fun/analyze_CIMP_all_CGIs_bis.R')
out <- analyze_CIMP_all_CGIs_bis(DiseaseList,CIMP.Number = 2)
```

**4)** We then analyze whether the methylation aberrations can be associated with transcriptomic or genetic variations: 


**4.A)** Can we assess CIMP from gene expression variations i.e CIMP=f(Gene Expression)?

We propose to tackle this problem using a sparse logistic regression with different formulations:

**i.** In the first case we predict the CIMP status for each tissue separately: 

```{r}
source('fun/predict_CIMP_GE_glmnet.R')
for (DiseaseName in DiseaseList)
{
        print(DiseaseName)
        out <- predict_CIMP_GE(DiseaseName, var.thresh=var.thresh, CIMP.Number=2, centered=T, scaled=T, intercept=T, n.folds=3, bootstrap=100, cores=10, log_exp=T, balanced=T)
}
```

**ii.** In the second case we compute a single classifier for all datasets:

```{r}
source('fun/predict_CIMP_GE_all.R')
out <- predict_CIMP_GE_all(DiseaseList, var.thresh=var.thresh, CIMP.Number=2, centered=T, scaled=T, intercept=T, n.folds=3, bootstrap=100, cores=10, log_exp=T, balanced=T)
```

**iii.** Finally, in the last case we relax the previous constraint (single classifier) by forcing each tissue-specific predictor to have the same non-zero coefficients but allowing the coefficients to vary:

```{r}
source('fun/predict_CIMP_GE_MT_par.R')
out <- predict_CIMP_GE_MT(DiseaseList, var.thresh=var.thresh, CIMP.Number=2, centered=T, scaled=T, intercept=T,  n.folds=3, bootstrap=100, cores=10, balanced=T)
```

**iv.** Summary of the results: 

```{r}
source('fun/analyze_predict_CIMP_GE_MT.R')
```

**4.B)** Analysis of the mutations associated with CIMP:

**i)** We analyzed the the association between CIMP and known reported mutations associated with tissue-specific CIMPs (_e.g_ _BRAF_, _KRAS_, _IDH1_, _IDH2_, _TET2_).

```{r}
source('fun/analyze_mutations.R')
Mutation.List <- c('BRAF','KRAS','IDH1','IDH2','TET2')

analyze_mutations(DiseaseList, Mutation.List=Mutation.List)
```

**ii)** We then also searched for other mutations significantly associated with CIMP in all diseases: 
```{r}
source('fun/analyze_mutations.R')
analyze_mutations(DiseaseList, Mutation.List=Mutation.List)
```

**5)** Survival analysis

```{r}
source('fun/compare_clinical.R')
for (DiseaseName in DiseaseList)
{
        print(DiseaseName)
        out <- compare_clinical(DiseaseName=DiseaseName, var.thresh=var.thresh, CIMP.Number= CIMP.Number)
}
```















