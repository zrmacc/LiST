# Purpose: Repeated Score Test
# Updated: 180413

#' Repeated Score Test for Normal Models
#' 
#' Individually tests \eqn{H_{0}:\beta_{j}=0} for each column of \code{X1}, 
#' adjusting for \code{X2}.
#' 
#' @param y Response vector.
#' @param X1 Numeric matrix of covariates whose regression coefficients are 
#'   constrained under the null.
#' @param X2 Numeric matrix of covariates whose regression coefficients are 
#'   unconstrained under the null.
#' @param tau Vaule of the scale parameter if known. If omitted, \eqn{\tau} is 
#'   estimated.
#' @return Vector of p-values, one for each column of \code{X1}, based on the 
#'   chi square distribution.
#'   
#' @importFrom plyr aaply
#' @importFrom stats model.matrix pchisq
#' @export

rScore.nlm = function(y,X1,X2,tau){
  # Missingness
  A = cbind(y,X2);
  aux = function(x){sum(is.na(x))>0};
  keep = !apply(A,MARGIN=1,FUN=aux);
  if(sum(!keep)>0){
    warning("Missing data detected in X2. These observations are excluded.")
    y = y[keep];
    X2 = X2[keep,];
  };
  # Fit null model
  if(missing(tau)){
    M0 = fitNorm(y=y,Z=X2,estT=T,t=1);
  } else {
    M0 = fitNorm(y=y,Z=X2,estT=F,t=tau);
  }
  # Extract residuals and scale residuals
  eT = M0$eT
  I22 = M0$Tau*M0$Ibb;
  # Function to calculate score statistics
  aux = function(x){
    I11 = sum(x^2);
    I12 = fastIP(A=x,B=X2);
    V = as.numeric(SchurC(I11=I11,I22=I22,I12=I12));
    a = as.numeric(fastIP(A=x,B=eT));
    Ts = a^2/(V*M0$Tau);
    return(Ts);
  }
  # Calculate score statistics
  S = aaply(.data=X1,.margins=2,.fun=aux);
  # Calculate p-value
  p = pchisq(q=S,df=1,lower.tail=F);
  return(p);
}