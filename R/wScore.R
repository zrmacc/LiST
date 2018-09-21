# Purpose: Weighted Score Test
# Updated: 180920

#' Weighted Score Test
#' 
#' Tests the hypothesis that a subset of the regression coefficients are fixed 
#' at a reference value, affording different weights to different coefficients. 
#' The weights are renormalized internally such that the observed test statistic
#' is one.
#' 
#' @param y Numeric response vector.
#' @param X Numeric model matrix.
#' @param L Logical vector, with as many entires as columns in the model matrix,
#'   indicating which columns have fixed coefficients under the null.
#' @param b10 Value of the regression coefficient for the selected columns under
#'   the null. Defaults to zero.
#' @param w Weights for the selected regression coefficients under the null.
#'   Defaults to the one vector.
#' @param method Either "asymptotic" or "perturbation".
#' @param B Score perturbations.
#' @param parallel Run in parallel? Must reigster parallel backend first.
#' 
#' @return A numeric vector containing the score statistic, the degrees of 
#'   freedom, and a p-value. 
#' 
#' @importFrom CompQuadForm davies
#' @importFrom foreach foreach '%dopar%' 
#' @importFrom stats rbinom
#' @export
#' 
#' @examples
#' \dontrun{
#' # See Vignette for Data Generation.
#' # Hypothesis test
#' L = c(F,F,F,T,T);
#' # Asymptotic p-value
#' wScore(y=y,X=X,L=L,w=c(2,1),method="asymptotic");
#' # Perturbation p-value
#' wScore(y=y,X=X,L=L,w=c(2,1),method="perturbation");
#' }

wScore = function(y,X,L,b10=NULL,w=NULL,method="asymptotic",B=1e3,parallel=F){
  ## Input checks
  if(!is.vector(y)){stop("A numeric vector is expected for y.")};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  if(!is.logical(L)){stop("A logical vector is expected for L.")};
  # Test specification
  q = ncol(X);
  k = sum(L);
  if(length(L)!=q){stop("L should have as many entries as columns in X.")};
  if(k==0){stop("At least 1 entry of L should be TRUE.")};
  if(k==q){stop("At least 1 entry of L should be FALSE.")};
  # Check for missingness
  Miss = sum(is.na(y))+sum(is.na(X));
  if(Miss>0){stop("Neither y nor Z, should contain missing values.")};
  # Null coefficient
  if(is.null(b10)){b10=rep(0,times=k)};
  if(length(b10)!=k){stop("b10 should contain a reference value for each TRUE element of L.")};
  # Weights
  if(is.null(w)){w=rep(1,times=k)};
  if(length(w)!=k){stop("w should contain a weight for each TRUE element of L.")};
  # Method
  Choices = c("asymptotic","perturbation");
  if(!(method%in%Choices)){stop("Select method from among asymptotic or perturbation.")};
  # Unparallelize
  if(!parallel){registerDoSEQ()};
  
  ## Partition target design
  Xa = X[,L,drop=F];
  Xb = X[,!L,drop=F];
  # Adjust response for fixed component
  y = y-as.numeric(MMP(Xa,b10));
  # Fit the null model
  M0 = fitOLS(y=y,X=Xb);
  # Calculate the score
  U = matIP(Xa,M0$Resid);
  # Score statistic, using raw weights
  Ts = as.numeric(matIP(A=U,B=(w^2)*U));
  # Renormalize weights
  w = w/sqrt(Ts);
  # Normalized score statistic
  Ts = 1;
  # Estimate p-value
  if(method=="asymptotic"){
    ## Asymptotic
    # Efficient information
    I11 = matIP(Xa,Xa);
    I12 = matIP(Xa,Xb);
    I22 = (M0$Ibb)*M0$V;
    V = SchurC(I11,I22,I12);
    # Xi
    Xi = M0$V*matQF(X=diag(w),A=V);
    # Eigenvalues
    lambda = eigen(x=Xi,symmetric=T,only.values=T)$values;
    lambda = round(lambda,digits=6);
    # Calculate p-value
    p = davies(q=Ts,lambda=lambda)$Qq;
  } else {
    # Perturbations
    R = foreach(i=1:(B-1),.combine=c) %dopar% {
      # Obs
      n = length(y);
      # Weights
      u = 2*rbinom(n=n,size=1,prob=0.5)-1;
      # Perturbed score
      Ub = matIP(Xa,u*M0$Resid);
      # Statistic
      Tb = as.numeric(matIP(A=Ub,B=(w^2)*Ub));
      return(Tb);
    }
    # Perturbation p
    p = (1+sum(R>=Ts))/(1+B);
  };
  # Output
  Out = c(Ts,k,p);
  names(Out) = c("Score","df","p");
  return(Out);
}