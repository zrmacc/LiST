#' @useDynLib NormalScoreTest
#' @importFrom Rcpp sourceCpp
NULL

#' Score Test for Normal Models
#' 
#' Performs the score test of the null hypothesis \eqn{\beta = 0} for a normal 
#' linear model of the form: \deqn{y = G\beta + X\alpha + \epsilon} Here 
#' \eqn{\epsilon} is distributed as \eqn{N(0,\tau R)}, where \eqn{\tau} is a 
#' variance component, and \eqn{R} is an \eqn{n x n} correlation structure.
#' 
#' @param G Matrix whose regression coefficient is zero in the null model.
#' @param X Matrix whose regression coefficient is unconstrained in the null 
#'   model.
#' @param y Response vector.
#' @param R Optional fixed correlation structure for variance component. Default
#'   is identity.
#' @return Vector including the score statistic, the degrees of freedom, and a 
#'   p-value based on the chi-square distribution.
#'   
#' @importFrom stats model.matrix pchisq
#' @importFrom matrixcalc is.diagonal.matrix
#' @export
#' 
#' @examples 
#' # Genotypes
#' g = NormalScoreTest::G[,1,drop=F];
#' # Phenotype
#' y = NormalScoreTest::Y[,1];
#' NormalScoreTest::NormScore(G=g,X=NormalScoreTest::X,y=y);

NormScore = function(G,X,y,R,tau.1=F){
  # Input formatting
  y = as.numeric(y);
  X = data.frame(X);
  G = data.frame(G);
  # Observations
  n = nrow(X);
  # Covariates
  k = ncol(G);
  q = ncol(X);
  # Check for correlation structure
  if(missing(R)){
    Ri = diag(n)
  } else if(matrixcalc::is.diagonal.matrix(R)){
    Ri = solve(R);
  } else {
    Ri = fastInv(R);
  }
  # Check dimensions
  check.dim = (nrow(G)==n)&(length(y)==n)&(nrow(Ri)==n);
  if(!check.dim){stop("Input dimensions are inconsistent")};
  # Design matrix for X
  Dx = model.matrix(~.,data=X);
  # Estimate alpha under H0
  alpha0 = Alpha0(X=Dx,Ri=Ri,y=y);
  # Estimate tau under H0
  tau0 = Tau0(X=Dx,Ri=Ri,y=y,a=alpha0);
  # Design matrix for G 
  Dg = model.matrix(~0+.,data=G);
  # Score statistic
  Ts = nScore(X=Dx,G=Dg,Ri=Ri,y=y,tau=tau0);
  # Calculate p-value
  p = pchisq(q=Ts,df=k,lower.tail=F);
  # Output
  Out = c(Ts,k,p);
  names(Out) = c("Score","df","p");
  return(Out);
}

#' Repeated Score Test
#' 
#' For each column of \eqn{G_{j}}, tests \eqn{H_{0}:\beta_{j}=0} in the model: 
#' \deqn{y = G_{j}\beta_{j}+X\alpha+\epsilon} See \code{\link{ScoreTest}} for 
#' model description.
#' 
#' @param G Matrix whose regression coefficients are zero under the null model.
#' @param X Matrix whose regression coefficient is unconstrained in the null 
#'   model.
#' @param y Response vector.
#' @param R Optional fixed correlation structure for variance component. Default
#'   is identity.
#' @param parallel Run in parallel? Must register parallel backend first.
#' @return Matrix with one row per column of G. Each row includes the score
#'   statistic, the degrees of freedom, and a p-value based on the chi-square
#'   distribution.
#' 
#' @importFrom foreach "%dopar%" foreach registerDoSEQ
#' @importFrom matrixcalc is.diagonal.matrix
#' @importFrom stats model.matrix pchisq
#' @export 
#' 
#' @examples 
#' # Genotype matrix
#' G = NormalScoreTest::G;
#' # Standardizing genotype matrix
#' Gs = scale(G);
#' # Phenotype
#' y = NormalScoreTest::Y[,1];
#' # Apply score test to each column of Gs
#' Results = NormalScoreTest::rNormScore(G=Gs,X=NormalScoreTest::X,y=y);

rNormScore = function(G,X,y,R,parallel=F){
  # Input formatting
  y = as.numeric(y);
  X = data.frame(X);
  G = data.frame(G);
  # Observations
  n = nrow(X);
  # Covariates
  k = ncol(G);
  q = ncol(X);
  # Check for correlation structure
  if(missing(R)){
    Ri = diag(n)
  } else if(matrixcalc::is.diagonal.matrix(R)){
    Ri = solve(R);
  } else {
    Ri = fastInv(R);
  }
  # Check dimensions
  check.dim = (nrow(G)==n)&(length(y)==n)&(nrow(Ri)==n);
  if(!check.dim){stop("Input dimensions are inconsistent")};
  # Estimate alpha under H0
  Dx = model.matrix(~.,data=X);
  alpha0 = Alpha0(X=Dx,Ri=Ri,y=y);
  # Estimate tau
  tau0 = Tau0(X=Dx,Ri=Ri,y=y,a=alpha0);
  if(!parallel){registerDoSEQ()};
  Out = foreach(i=1:k,.combine=rbind,.inorder=T) %dopar%{
    # Extract covariate
    g = G[,i];
    # Design matrix for G
    Dg = model.matrix(~0+.,data=data.frame(g));
    # Score statistic
    Ts = nScore(X=Dx,G=Dg,Ri=Ri,y=y,tau=tau0);
    # Calculate p-value
    p = pchisq(q=Ts,df=1,lower.tail=F);
    # Results
    Res = c(Ts,1,p);
    names(Res) = c("Score","df","p");
    return(Res);
  }
  return(Out)
}

#' Kernelized Score Test for Linear Models
#' 
#' Performs the score test of the null hypothesis \eqn{\beta = 0} for a normal 
#' linear model of the form: \deqn{y = G\beta + X\alpha + \epsilon} Here 
#' \eqn{\epsilon} is distributed as \eqn{N(0,\tau R)}, where \eqn{\tau} is a 
#' variance component, and \eqn{R} is an \eqn{n x n} correlation structure. The 
#' kernelized score statistic takes the form \deqn{T_{K} =
#' \tilde{\epsilon}'K\tilde{\epsilon}} Here \eqn{\tilde{\epsilon}} is the
#' residual estimated under \eqn{H_{0}}, and \eqn{K = WG(WG)'} is the kernel matrix.
#' 
#' @param G Matrix whose regression coefficient is zero in the null model.
#' @param W Weight matrix.
#' @param X Matrix whose regression coefficient is unconstrained in the null 
#'   model.
#' @param y Response vector.
#' @param R Optional fixed correlation structure for variance component. Default
#' is identity.
#' @return List including the \code{Score}, the \code{Eigenvalues} of the
#'   mixutre distribution, and the \code{p} value.
#' 
#' @importFrom CompQuadForm davies
#' @importFrom matrixcalc is.diagonal.matrix is.positive.semi.definite
#' @export
#' 
#' @examples 
#' # Genotypes
#' G = NormalScoreTest::G[,1:10];
#' maf = apply(G,MARGIN=2,FUN=mean)/2;
#' Gs = scale(G);
#' # Weights
#' W = round(diag((dbeta(x=maf,1,25))),digits=6);
#' # Phenotype
#' y = NormalScoreTest::Y[,1];
#' Results = NormalScoreTest::kNormScore(G=Gs,W=W,X=NormalScoreTest::X,y=y);

kNormScore = function(G,W,X,y,R){
  # Input Checks
  k = ncol(G);
  if(nrow(W)!=k){stop("Dimension of weight matrix should match ncol(G).")}
  # Kernel matrix
  if(! matrixcalc::is.diagonal.matrix(W)){stop("Diagonal weight matrix is expected.")};
  L = G %*% W;
  # Input formatting
  y = as.numeric(y);
  X = data.frame(X);
  G = data.frame(G);
  # Observations
  n = nrow(X);
  # Covariates
  q = ncol(X);
  # Check for correlation structure
  if(missing(R)){
    R = Ri = diag(n);
  } else if(matrixcalc::is.diagonal.matrix(R)){
    Ri = solve(R);
  } else {
    Ri = fastInv(R);
  }
  # Check dimensions
  check.dim = (nrow(G)==n)&(length(y)==n)&(nrow(Ri)==n);
  if(!check.dim){stop("Input dimensions are inconsistent")};
  # Number of covariates
  q = ncol(X);
  # Estimate alpha under H0
  Dx = model.matrix(~.,data=X);
  alpha0 = Alpha0(X=Dx,Ri=Ri,y=y);
  # Estimate tau
  tau0 = Tau0(X=Dx,Ri=Ri,y=y,a=alpha0);
  # Score statistic
  S = knScore(X=Dx,L=L,Ri=Ri,y=y,tau=tau0);
  Tk = S$Tk;
  Q = S$Q;
  # Eigenvalues
  U = matIP(X=L,A=Q);
  lambda = eigen(x=U,symmetric=T,only.values=T)$values;
  lambda = round(lambda,digits=6);
  # Calculate p-value
  p = CompQuadForm::davies(q=Tk,lambda=lambda)$Qq;
  # Output
  Out = list("Score"=Tk,"Eigenvalues"=lambda,"p"=p);
  return(Out);
}