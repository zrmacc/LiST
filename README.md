---
title: "README"
author: "Zachary McCaw"
date: "2018-09-21"
output: 
  html_document: 
    keep_md: TRUE
--- 

# Package Vignette


# Contents

* [Setting](#setting)
* [Example Data](#example-data)
* [Score Tests](#score-test)

# Setting
Consider the standard linear model:

$$
y = X\beta + \epsilon,\ \epsilon\sim(0,\sigma^2)
$$

Partition the regression coefficient $\beta=(\beta_{1},\beta_{2})$. Suppose interest lies in testing whether $\beta_{1}$ is fixed at a reference value. This package implements score tests of $H_{0}:\beta_{1}=\beta_{1}^{\dagger}$. 

# Example Data

Below, data are simulated for $n=10^{3}$ subjects. The outcome $y_{i}\in\mathbb{R}$ is normally distributed with unit variance. The subject-specific means depend on an intercept and four covariates. The function `Fit.LinReg` provides an implementation of linear regression. 


```r
library(LiST);
# Generate data
set.seed(102);
# Observations
n = 1e3;
# Design
X = cbind(1,matrix(rnorm(4*n),nrow=n));
colnames(X) = c("int",paste0("x",seq(1:4)));
# Regression coefficient
b = c(1,-0.1,0.2,-0.1,0);
# Linear predictor
h = as.numeric(X%*%b);
# Outcome
y = h+rnorm(n);
# Fitted linear model
M = Fit.LinReg(y=y,X=X);
show(M);
```

```
## Fitted Linear Model
## Estimated Coefficients:
##   Coeff   Point     SE       L       U         p
## 1   int  1.0200 0.0315  0.9600  1.0800 2.51e-231
## 2    x1 -0.1430 0.0328 -0.2070 -0.0786  1.34e-05
## 3    x2  0.1780 0.0310  0.1170  0.2380  1.01e-08
## 4    x3 -0.0839 0.0317 -0.1460 -0.0218  8.09e-03
## 5    x4 -0.0130 0.0322 -0.0761  0.0500  6.85e-01
```

# Score Tests

## Standard Test

`Score` conducts a standard score test, using either the asymptotic variance or perturbation to estimate a $p$-value. The hypothesis test is specified using a logical vector `L`. Elements of `L` that are constrained under the null are set to `TRUE`, those which are estimated under the null are set to `FALSE`. At least one element of `L` must be `TRUE` (a test must be specified), and at least one element of `L` must be `FALSE` (the null model must be estimable). Below, hypotheses are tested on the example data.

* The first test assesses $H_{0}:\beta_{4}=0$, which is true.
* The second test assesses $H_{0}:\beta_{2}=0$, which is false.
* The third test assesses $H_{0}:\beta_{1}=\beta_{3}=-0.1$, which is true. 


```r
cat("Score test of b4=0\n");
Score(y=y,X=X,L=c(F,F,F,F,T),method="asymptotic");
Score(y=y,X=X,L=c(F,F,F,F,T),method="perturbation");
cat("\n");
cat("Score test of b2=0\n");
Score(y=y,X=X,L=c(F,F,T,F,F),method="asymptotic");
Score(y=y,X=X,L=c(F,F,T,F,F),method="perturbation");
cat("\n");
cat("Score test of b1=b3=-0.1\n");
Score(y=y,X=X,L=c(F,T,F,T,F),b10=c(-0.1,-0.1),method="asymptotic");
Score(y=y,X=X,L=c(F,T,F,T,F),b10=c(-0.1,-0.1),method="perturbation");
cat("\n");
```

```
## Score test of b4=0
##     Score        df         p 
## 0.1647543 1.0000000 0.6848166 
##     Score        df         p 
## 0.1647543 1.0000000 0.6543457 
## 
## Score test of b2=0
##        Score           df            p 
## 3.180192e+01 1.000000e+00 1.707236e-08 
##        Score           df            p 
## 31.801923115  1.000000000  0.000999001 
## 
## Score test of b1=b3=-0.1
##     Score        df         p 
## 2.0582910 2.0000000 0.3573121 
##     Score        df         p 
## 2.0582910 2.0000000 0.3396603
```

## Weighted Score Test

`wScore` conducts a score test that affords different components of the hypothesis different weights. The following assesses $H_{0}:\beta_{3}=\beta_{4}=0$, equivalently $H_{0}:(\beta_{3}=0)\land(\beta_{4}=0)$, giving different weights to the different components. The first component $(\beta_{3}=0)$ is false, while the second component $(\beta_{4}=0)$ is true. The overall hypothesis is false. Whether or not the weighted test rejects depends on which component of the hypothesis is afforded more weight. 


```r
cat("Weighted score test of b3=b4=0\n");
wScore(y=y,X=X,L=c(F,F,F,T,T),w=c(10,1),method="asymptotic");
wScore(y=y,X=X,L=c(F,F,F,T,T),w=c(10,1),method="perturbation");
cat("\n");
cat("Weighted score test of b3=b4=0\n");
wScore(y=y,X=X,L=c(F,F,T,F,T),w=c(1,10),method="asymptotic");
wScore(y=y,X=X,L=c(F,F,T,F,T),w=c(1,10),method="perturbation");
```

```
## Weighted score test of b3=b4=0
##       Score          df           p 
## 1.000000000 2.000000000 0.008175889 
##       Score          df           p 
## 1.000000000 2.000000000 0.006993007 
## 
## Weighted score test of b3=b4=0
##     Score        df         p 
## 1.0000000 2.0000000 0.4700908 
##     Score        df         p 
## 1.0000000 2.0000000 0.4445554
```
