# Purpose: Repeated Score Test
# Updated: 180409

#' Repeated Score Test for Normal Models
#' 
#' Individually tests \eqn{H_{0}\beta_{j}=0} for each column of \code{X1}, 
#' adjusting for \code{X2}.
#' 
#' @param y Response vector.
#' @param X1 Numeric matrix of covariates whose regression coefficients are 
#'   constrained under the null.
#' @param X2 Numeric matrix of covariates whose regression coefficients are 
#'   unconstrained under the null.
#' @param tau Vaule of the scale parameter if known. If omitted, \eqn{\tau} is 
#'   estimated.
#' @param K Optional fixed correlation structure for variance component. Default
#'   is identity.
#' @return Vector of p-values, one for each column of \code{X1}, based on the 
#'   chi square distribution.
#'   
#' @importFrom plyr aaply
#' @importFrom stats model.matrix pchisq
#' @export

rScore.nlm = function(y,X1,X2,tau,K){
  # Missingness
  A = cbind(y,X2);
  aux = function(x){sum(is.na(x))>0};
  keep = !apply(A,MARGIN=1,FUN=aux);
  if(sum(!keep)>0){
    warning("Missing data detected in X2. These observations are excluded.")
    y = y[keep];
    X2 = X2[keep,];
  };
  # Calculate score statistic
  if(missing(tau)){
    # Note: when tau is estimate, Ts is independent of K
    aux = function(x){normScore(y=y,X1=x,b=0,X2=X2,estT=T,t=1,useK=F,K=1)};
  } else {
    if(missing(K)){
      aux = function(x){normScore(y=y,X1=x,b=0,X2=X2,estT=F,t=tau,useK=F,K=1)};
    } else {
      aux = function(x){normScore(y=y,X1=X1,b=0,X2=X2,estT=F,t=tau,useK=T,K=K)};
    }
  }
  # Calculate score statistics
  S = aaply(.data=X1,.margins=2,.fun=aux);
  # Calculate p-value
  p = pchisq(q=S,df=1,lower.tail=F);
  return(p);
}