#' Simulated Genotypes
#' 
#' Simulated genotypes for 100 subjects at 1000 loci. Genotypes were simulated 
#' with use of hapgen2, and are based on the haplotype structure of human 
#' chromosome one in the CEU population of the 1000 Genomes Project.
#' @format An integer matrix with 1000 rows and 100 columns
#' \describe{ 
#'    \item{s1-s100}{s[i] is an integer vector of minor allele counts for the ith subject.}
#' }
"G"

#' Simulated Covariates
#' 
#' Age, sex, and two principal components of the centered and scaled subject by
#' locus genotype matrix. Age was drawn from a gamma distribution with mean 50
#' and variance 10. Sex was drawn from a Bernoulli distribuiton with expectation
#' 1/2.
#' @format A numeric matrix with 100 rows and 4 columns
#' \describe{
#'    \item{Age}{Age.}
#'    \item{Sex}{Sex.}
#'    \item{PC1}{First principal component.}
#'    \item{PC2}{First principal component.}
#' }
"X"

#' Simulated Phenotype
#' 
#' Phenotypes simulated under the null hypothesis of no
#' genotypic effect. A subject specific linear predictor was calculated based on
#' age, sex, pc1, and pc2. Residuals were drawn from a N(0,1) distribution. 
#' @format A numeric matrix with 100 rows and 1 columns
#' \describe{
#'    \item{YN}{Normal phenotype.}
#' }
"Y"

