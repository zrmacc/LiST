# Purpose: Score tests for normal linear models
# Updated: 180408

#' @useDynLib NST
#' @importFrom Rcpp sourceCpp
NULL

#' Score Test for Normal Models
#' 
#' Performs a score test of \eqn{H_{0}:\beta_{1}=\beta_{10}}.
#' 
#' @param y Response vector.
#' @param Z Numeric model matrix.
#' @param L Logical vector, with as many entires as columns in the model matrix,
#'   indicating which columns have fixed coefficients under the null.
#' @param b0 Value of the regression coefficient for the selected columns under 
#'   the null. Defaults to zero.
#' @param tau Vaule of the scale parameter if known. If omitted, \eqn{\tau} is 
#'   estimated.
#' @param K Optional fixed correlation structure for variance component. Default
#'   is identity.
#' @return Vector including the score statistic, the degrees of freedom, and a 
#'   p-value based on the chi squared distribution.
#'   
#' @importFrom stats model.matrix pchisq
#' @importFrom matrixcalc is.diagonal.matrix
#' @export

Score.nlm = function(y,Z,L,b0,K,tau){
  # Input checks
  if(!is.logical(L)){stop("L should be a logical vector.")};
  if(length(L)!=ncol(Z)){stop("L should have as many entries as columns in Z.")};
  if(sum(L)==0){stop("At least 1 entry of L should be TRUE.")}
  if(sum(L)==length(L)){stop("At least 1 entry of L should be FALSE.")};
  # Missingness
  A = cbind(y,Z);
  aux = function(x){sum(is.na(x))>0};
  keep = !apply(A,MARGIN=1,FUN=aux);
  if(sum(!keep)>0){
    warning("Missing data detected. These observations are excluded.")
    y.t = y.t[keep];
    y.s = y.s[keep];
    Z.t = Z.t[keep,];
    Z.s = Z.s[keep,];
  };
  # Null coefficient
  if(missing(b0)){b0=rep(0,times=sum(L))};
  # Partition model matrix
  X1 = Z[,L,drop=F];
  k = ncol(X1);
  X2 = Z[,!L,drop=F];
  # Calculate score statistic
  if(missing(tau)){
    # Note: when tau is estimate, Ts is independent of K
    Ts = normScore(y=y,X1=X1,b=b0,X2=X2,estT=T,t=1,useK=F,K=1);
  } else {
    if(missing(K)){
      Ts = normScore(y=y,X1=X1,b=b0,X2=X2,estT=F,t=tau,useK=F,K=1);
    } else {
      Ts = normScore(y=y,X1=X1,b=b0,X2=X2,estT=F,t=tau,useK=T,K=K);
    }
  }
  # Calculate p-value
  p = pchisq(q=Ts,df=k,lower.tail=F);
  # Output
  Out = c(Ts,k,p);
  names(Out) = c("Score","df","p");
  return(Out);
}

#' Weighted Score Test for Normal Models
#' 
#' Performs a weighted score test of \eqn{H_{0}:\beta_{1}=\beta_{10}} using
#' the test statistic \deqn{T_{W}=S'WWS} Here \eqn{S} is the score vector and
#' \eqn{W} is a summetric weight matrix.
#' 
#' @param y Response vector.
#' @param Z Numeric model matrix.
#' @param W Numeric weight matrix.
#' @param L Logical vector, with as many entires as columns in the model matrix,
#'   indicating which columns have fixed coefficients under the null.
#' @param b0 Value of the regression coefficient for the selected columns under 
#'   the null. Defaults to zero.
#' @param tau Value of the scale parameter if known. If omitted, \eqn{\tau} is 
#'   estimated.
#' @param K Optional fixed correlation structure for variance component. Default
#'   is identity.
#' @return Vector including the score statistic, the degrees of freedom, and a 
#'   p-value based on the chi squared distribution.
#' 
#' @importFrom CompQuadForm davies
#' @importFrom matrixcalc is.diagonal.matrix is.positive.semi.definite
#' @export

wScore.nlm = function(y,Z,W,L,b0,tau,K){
  # Input checks
  if(!is.logical(L)){stop("L should be a logical vector.")};
  if(length(L)!=ncol(Z)){stop("L should have as many entries as columns in Z.")};
  if(sum(L)==0){stop("At least 1 entry of L should be TRUE.")}
  if(sum(L)==length(L)){stop("At least 1 entry of L should be FALSE.")};
  if(ncol(W)!=length(y)){stop("W should have the same dimensions as the number of observations.")}
  # Missingness
  A = cbind(y,Z);
  aux = function(x){sum(is.na(x))>0};
  keep = !apply(A,MARGIN=1,FUN=aux);
  if(sum(!keep)>0){
    warning("Missing data detected. These observations are excluded.")
    y.t = y.t[keep];
    y.s = y.s[keep];
    Z.t = Z.t[keep,];
    Z.s = Z.s[keep,];
  };
  # Null coefficient
  if(missing(b0)){b0=rep(0,times=sum(L))};
  # Partition model matrix
  X1 = Z[,L,drop=F];
  k = ncol(X1);
  X2 = Z[,!L,drop=F];
  # Calculate score statistic
  if(missing(tau)){
    # Note: when tau is estimate, Ts is independent of K
    S = wNormScore(y=y,X1=X1,b=b0,X2=X2,estT=T,t=1,useK=F,K=1,W=W);
  } else {
    if(missing(K)){
      S = wNormScore(y=y,X1=X1,b=b0,X2=X2,estT=F,t=tau,useK=F,K=1,W=W);
    } else {
      S = wNormScore(y=y,X1=X1,b=b0,X2=X2,estT=F,t=tau,useK=T,K=K,W=W);
    }
  }
  # Eigenvalues
  Xi = matQF(X=W,A=S$Q);
  lambda = eigen(x=Xi,symmetric=T,only.values=T)$values;
  lambda = round(lambda,digits=6);
  # Calculate p-value
  p = CompQuadForm::davies(q=S$Tw,lambda=lambda)$Qq;
  # Output
  Out = list("Score"=S$Tw,"Eigenvalues"=lambda,"p"=p);
  return(Out);
}