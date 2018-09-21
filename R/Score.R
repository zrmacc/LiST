# Purpose: Standard score test
# Updated: 180920

#' Standard Score Test
#' 
#' Tests the hypothesis that a subset of the regression coefficients are fixed 
#' at a reference value. Specifically, let \eqn{\beta} denote the regression 
#' coefficient. Partition \eqn{\beta=(\beta_{1},\beta_{2})}. Suppose that 
#' interest lies in testing that \eqn{\beta_{1}} is fixed at \eqn{\beta_{10}}. 
#' \code{Score} performs a score test of \eqn{H_{0}:\beta_{1}=\beta_{10}}. The
#' test is specified using a logical vector \code{L}, with as many entries as
#' columns in the model matrix \code{X}. The values of \code{L} set to \code{T}
#' are constrained under the null, while values of \code{L} set to \code{F} are
#' estimated under the null.
#' 
#' @param y Numeric response vector.
#' @param X Numeric model matrix.
#' @param L Logical vector, with as many entires as columns in the model matrix,
#'   indicating which columns have fixed coefficients under the null.
#' @param b10 Value of the regression coefficient for the selected columns under
#'   the null. Defaults to zero.
#' @param method Either "asymptotic" or "perturbation".
#' @param B Score perturbations.
#' @param parallel Run in parallel? Must reigster parallel backend first. 
#' 
#' @return A numeric vector containing the score statistic, the degrees of 
#'   freedom, and a p-value. 
#' 
#' @importFrom foreach registerDoSEQ foreach '%dopar%'
#' @importFrom stats rbinom pchisq
#' @export
#' 
#' @examples
#' \dontrun{
#' # See Vignette for Data Generation.
#' # Hypothesis test
#' L = c(F,F,F,F,T);
#' # Asymptotic p-value
#' Score(y=y,X=X,L=L,method="asymptotic");
#' # Perturbation p-value
#' Score(y=y,X=X,L=L,method="perturbation");
#' }

Score = function(y,X,L,b10=NULL,method="asymptotic",B=1e3,parallel=F){
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
  # Effincient information
  I11 = matIP(Xa,Xa);
  I12 = matIP(Xa,Xb);
  I22 = (M0$Ibb)*M0$V;
  E = SchurC(I11,I22,I12);
  # Score statistic
  Ts = (1/M0$V)*as.numeric(matQF(U,matInv(E)));
  # Estimate p-value
  if(method=="asymptotic"){
    # Asymptotic
    p = pchisq(q=Ts,df=k,lower.tail=F);
  } else {
    # Perturbations
    R = foreach(i=1:(B-1),.combine=c) %dopar% {
      # Obs
      n = length(y);
      # Weights
      w = 2*rbinom(n=n,size=1,prob=0.5)-1;
      # Perturbed score
      Ub = matIP(Xa,w*M0$Resid);
      # Statistic
      Tb = (1/M0$V)*as.numeric(matQF(Ub,matInv(E)));
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