README
================
Zachary McCaw
11/07/2017

Package Vignette
================

Contents
--------

-   [Standard Score Test](#standard-score-test)

Standard Score Test
-------------------

Consider the linear model:

*y* = *X**α* + *G**β* + *ϵ*,   *ϵ* ∼ *N*(0, *τ**R*)
 Here *β* is the regression parameter of interest, *α* is a nuisance regression parameter, *τ* is a variance component, and *R* is a fixed correlation structure. The function `lmScoreTest` provides a score test of *H*<sub>0</sub> : *β* = 0. In particular, the score statistic is calculated as:

$$
T\_{S} = S\_{\\beta}(\\beta=0,\\alpha=\\tilde{\\alpha})'I\_{\\beta\\beta'|\\alpha}^{-1}S\_{\\beta}(\\beta=0,\\alpha=\\tilde{\\alpha})'
$$
 Here $\\tilde{\\alpha}$ solves the equation *S*<sub>*α*</sub>(*β* = 0, *α*)=0. A *p*-value for *T*<sub>*S*</sub> is estimated by comparison with the *χ*<sub>dim(*β*)</sub><sup>2</sup> distribution.

``` r
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

    ## Estimated size at alpha=0.05:
    ## [1] 0.037
