# README
Zachary McCaw  
11/07/2017  

# Package Vignette




## Contents

* [Standard Score Test](#standard-score-test)

## Standard Score Test
Consider the linear model:

$$y=X\alpha+G\beta+\epsilon,\;\;\epsilon\sim N(0,\tau R)$$

Here $\beta$ is the regression parameter of interest, $\alpha$ is a nuisance regression parameter, $\tau$ is a variance component, and $R$ is a fixed correlation structure. The function `lmScoreTest` provides a score test of $H_{0}:\beta = 0$. In particular, the score statistic is calculated as:

$$T_{S}=S_{\beta}(\beta=0,\alpha=\tilde{\alpha})'I_{\beta\beta'|\alpha}^{-1}S_{\beta}(\beta=0,\alpha=\tilde{\alpha})'$$

Here $\tilde{\alpha}$ solves the equation $S_{\alpha}(\beta=0,\alpha) = 0$. A $p$-value for $T_{S}$ is estimated by comparison with the $\chi_{\dim(\beta)}^{2}$ distribution. 


```r
# Phenotype vector
Y = ScoreTest::Y;
# Covariate matrix
X = ScoreTest::X;
# Genotype matrix
G = ScoreTest::G;
# Standardizing genotype matrix
Gs = scale(t(G));
# Function to apply test across the loci in G 
aux = function(g){
  S = ScoreTest::lmScoreTest(G=g,X=X,y=Y);
  return(S["p"])
  }
P = apply(Gs,MARGIN=2,FUN=aux);
# Estimated size
cat("Estimated size at alpha=0.05:\n");
mean(P<=0.05);
```

```
## Estimated size at alpha=0.05:
## [1] 0.037
```
