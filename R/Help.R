# Purpose: Package documentation
# Updated: 180919

#' @useDynLib LiST
#' @importFrom Rcpp sourceCpp
NULL

#' LiST: Linear Score Tests
#' 
#' Implementation of score tests for linear models. \code{\link{Score}} performs a basic
#' score test, providing either an asymptotic or perturbation p-value. \code{\link{wScore}}
#' performs a weighted score test. \code{\link{rScore}} performs repeated association tests.
#' \code{\link{pScore}} performs score tests on a model previously fit using \code{\link{Fit.LinReg}}.
#' 
#' @author Zachary R. McCaw
#' @docType package
#' @name LiST-help
NULL