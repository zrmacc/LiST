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
  # Fit null model
  if(missing(tau)){
    # Note: when tau is estimate, Ts is independent of K
    M0 = fitNorm(y=y,Z=Z,estT=T,t=1,useK=F,K=1);
  } else {
    if(missing(K)){
      M0 = fitNorm(y=y,Z=Z,estT=F,t=tau,useK=F,K=1);
    } else {
      M0 = fitNorm(y=y,Z=Z,estT=F,t=tau,useK=T,K=K);
    }
  }
  # Extract and scale residuals
  eT = (M0$eT/M0$Tau);
  Q = M0$Q;
  # Function to calculate score statistics
  aux = function(x){
    a = as.numeric(fastIP(A=x,B=eT));
    b = vecQF(x=x,A=Q);
    Ts = a^2/b;
    return(Ts);
  }
  # Calculate score statistics
  S = aaply(.data=X1,.margins=2,.fun=aux);
  # Calculate p-value
  p = pchisq(q=S,df=1,lower.tail=F);
  return(p);
}