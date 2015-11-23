---
layout: teaching
authors: Matahi MOARII, Marine Le MORVAN.
title: bio-info, SVM and Graph-kernels
location: Mines Paristech
categories: teaching-2015
---

# Introduction

This page is the practical session of the "Support Vector Machines" module taught by [Chloé-Agathe Azencott](http://cazencott.info).

# Material

The practical session is done using [R](http://cran.r-project.org). For a quick tutorial, follow this [link](http://pages.pomona.edu/~jsh04747/courses/R.pdf).

# Practical Session

## Necessary support:

- The lecture [slides](http://cazencott.info/dotclear/public/lectures/2015-S1234-svm.pdf)
- The practical session R and the datasets [file](http://matahi.github.io/downloads/TPMines2015.zip) 

### Part A: Linear dataset.

**A.1/A.2)** Visualization of the first dataset:

```{r}
load('linear1.RData')
```

3 new variables are in the environment now (type `ls()` to check):

- linear1.train : a matrix of size n1x3 (first column is *x1*, second column is *x2* and third column is the output *y*). 
- linear1.test.input : a matrix of size n2 x 2 (same as the train but without the output *y*).
- linear1.data : the combined dataset (for visualization).

Then visualize the dataset:

```{r}
require( 'ggplot2' )
qplot( data=linear1.data, x.1, x.2, colour=factor(y), 
				    shape=factor(train) )
```

**Question 1:** How can you characterise this dataset ?

Given what we learned during the lecture, we train a **linear SVM** on the training set which will be our predictor for the testing set.


**A.3/A.4/A.5)** Training the SVM:

```{r}
### A.3) Train a linear SVM
require( 'kernlab' )
linear1.svm <- ksvm( y ~ ., data=linear1.train, type='C-svc', 
						kernel='vanilladot',
						C=100, scale=c() )

### A.4) Plot the model
plot( linear1.svm, data=linear1.train )

### A.5) Adding points of test on the graph
points( linear1.test.input[ sample.int(nrow(linear1.test.input),10), ], pch=4 )
```

**Question 2:** What are the black points in the figure ?

**A.6/A.7** Testing on another dataset:


```{r}
### A.6) Prediction
linear1.prediction <- predict( linear1.svm, linear1.test.input )

### A.7) Look at accuracy
load('linear1Sol.RData')
# contains linear1.test.output

print(paste0('Accuracy: ', floor(100*sum( linear1.prediction == linear1.test.output )/length(linear1.test.output)), '%'))
```

Let's consider a dataset a bit more complex:

**A.9 to A.14)** This is the exact same thing as the first part except that the dataset is non-separable.

**Question 3:** Is the accuracy a sufficient method to assess the performance of a model?

Let's discuss a few ways to improve the assessment of the performance:

**A.15)** Separate positive and negative examples :

```{r}
### A.15) A confusion matrix gives more information than just accuracy
print('Confusion Matrix: ');print(table( linear2.prediction, linear2.test.output, dnn= c("prediction","reality") ))
```

**A.16)** ROC Curves (See [wikipedia](http://en.wikipedia.org/wiki/Receiver_operating_characteristic)):

```{r}
linear2.prediction.score <- predict( linear2.svm, linear2.test.input, type='decision' )

require( 'ROCR' )

## ROC
linear2.roc.curve <- performance( prediction( linear2.prediction.score, linear2.test.output ), measure='tpr', x.measure='fpr' )
plot( linear2.roc.curve )
```

### Part B: Non-Linear dataset.

**B.1/B.2/B.3)** Trying linear SVM on a dataset where it is not appropriate.

Let's try a different kernel

**B.4)** RBF Kernel:

```{r}
nonlinear.svm <- ksvm( y ~ ., data=nonlinear.train, type='C-svc', 
			      kernel='rbf', kpar=list(sigma=1), 
			      C=100, scale=c() )
plot( nonlinear.svm, data=nonlinear.train )
```

**Question 4:** Recall what is the parameter `C` in the svm. In the following we will see what happens when we vary `C`.

**B.5)** Impact of `C`

```{r}
require('manipulate')
manipulate( plot( ksvm( y ~ ., data=nonlinear.train, type='C-svc', 
			kernel=k, C=2^c.exponent, scale=c() ),
 			data=nonlinear.train ), 
		   c.exponent=slider(-10,10),
                   k=picker('Gaussian'='rbfdot', 'Linear'='vanilladot', 
		            'Hyperbolic'='tanhdot','Spline'='splinedot', 
			    'Laplacian'='laplacedot') )
```

**B.6)** Visualization of the impact of `C` on the prediction accuracy (**Bias-Variance Tradeoff**):

```{r}
### B.6) Bias-Variance Tradeoff
BiasVarianceTradeoff <- function( dataset, cross=10, c.seq=2^seq(-10, 10), ... ) {
  err <- sapply( c.seq, function( c ) 
                {
                        cross( ksvm( y ~ ., data=dataset, C=c, cross=cross, ...) )
                })
  return(data.frame( c=c.seq, error=err ))
}

qplot( c, error, data=BiasVarianceTradeoff( nonlinear.train, type='C-svc', kernel='rbfdot' ), geom='line', log='x' )
```

**Question 5:** How to choose `C` ?


### Part C: Acute lymphoblastic leukemia dataset.

In this part, we work on a real public dataset of *Acute lymphoblastic leukemia (ALL)* patients.

**C.1)** We are interested in classifying leukemia patients into two classes: *B-cell ALL* vs *T-cell ALL*  because this classification can have an impact on the **patient's prognosis** or **its response to a given treatment**.

**C.2)** When we add supplementary clinical information, we have more than 2 classes (*B1*, *B2*,...).  

**Question 6:** How can you derive from what you learned a SVM classifier that can predict more than 2 classes?

