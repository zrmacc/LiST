#' @useDynLib ScoreTest
#' @importFrom Rcpp sourceCpp
NULL

#' Score Test for Linear Models
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
#' @param tau.1 Estimate \eqn{\tau} using the full model? If false, \eqn{\tau}
#'   is estimated under \eqn{H_{0}}.
#' @return Vector including the score statistic, the degrees of freedom, and a 
#'   p-value based on the chi-square distribution.
#'   
#' @importFrom stats model.matrix pchisq
#' @importFrom matrixcalc is.diagonal.matrix
#' @export
#' 
#' @examples 
#' ScoreTest::ScoreTest(G=ScoreTest::G[1,],X=ScoreTest::X,y=ScoreTest::Y);

ScoreTest = function(G,X,y,R,tau.1=F){
  # Input formatting
  X = data.frame(X);
  G = data.frame(G);
  # Observations
  n = nrow(X);
  # Number of covariates
  k = ncol(G);
  q = ncol(X);
  # Check for correlation structure
  if(missing(R)){
    Ri = diag(n)
  } else if(is.diagonal.matrix(R)){
    Ri = solve(R);
  } else {
    Ri = fastInv(R);
  }
  # Estimate alpha under H0
  Dx = model.matrix(~.,data=X);
  alpha0 = Alpha0(X=Dx,Ri=Ri,y=y);
  # Null residuals
  e0 = (y - Dx %*% alpha0);
  # Estimate tau
  if(tau.1){
    D = model.matrix(~.,data=data.frame(X,G))
    tau = Tau1(D=D,Ri=Ri,y=y);
  } else {
    tau = sum(e0^2)/(n-q);
  }
  # Design matrix for G
  Dg = model.matrix(~0+.,data=G);
  # Information matrix for (beta,alpha)
  J = Info(X=Dx,G=Dg,Ri=Ri,tau=tau);
  # Efficient information for beta
  K = SchurC(Igg=J$Ibb,Ihh=J$Iaa,Igh=J$Iba,inv=T);
  # Score for beta under H0
  S = ScoreB(e0=e0,G=Dg,Ri=Ri,tau=tau);
  # Score statistic
  Ts = qForm(K=K,s=S);
  # Calculate p-value
  p = pchisq(q=Ts,df=k,lower.tail=F);
  # Output
  Out = c(Ts,k,p);
  names(Out) = c("Score","df","p");
  return(Out);
}

#' Serial Score Test
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
#' @importFrom stats model.matrix pchisq
#' @importFrom foreach "%dopar%" foreach registerDoSEQ
#' @export 
#' 
#' @examples 
#' # Genotype matrix
#' G = ScoreTest::G;
#' # Standardizing genotype matrix
#' Gs = scale(t(G));
#' # Apply score test to each column of Gs
#' Results = ScoreTest::sScoreTest(G=Gs,X=ScoreTest::X,y=ScoreTest::Y);

sScoreTest = function(G,X,y,R,parallel=F){
  # Input formatting
  X = data.frame(X);
  G = data.frame(G);
  # Observations
  n = nrow(X);
  # Number of covariates
  k = ncol(G);
  q = ncol(X);
  # Check for correlation structure
  if(missing(R)){
    Ri = diag(n)
  } else if(is.diagonal.matrix(R)){
    Ri = solve(R);
  } else {
    Ri = fastInv(R);
  }
  # Estimate alpha under H0
  Dx = model.matrix(~.,data=X);
  alpha0 = Alpha0(X=Dx,Ri=Ri,y=y);
  # Null residuals
  e0 = (y - Dx %*% alpha0);
  # Estimate tau
  tau = sum(e0^2)/(n-q);
  # Information for alpha
  Iaa = (1/tau) * matIP(X=Dx,A=Ri);
  # Apply score test to columns of G
  if(!parallel){registerDoSEQ()};
  i = NULL;
  Out = foreach(i=1:k,.combine=rbind,.inorder=T) %dopar%{
    # Extract covariate
    g = G[,i];
    # Design matrix for G
    Dg = model.matrix(~0+.,data=data.frame(g));
    # Information for beta
    Ibb = (1/tau) * matIP(X=Dg,A=Ri);
    Iba = (1/tau) * t(Dg) %*% data.matrix(Dx);
    # Efficient information for beta
    K = SchurC(Igg=Ibb,Ihh=Iaa,Igh=Iba,inv=T);
    # Score for beta under H0
    S = ScoreB(e0=e0,G=Dg,Ri=Ri,tau=tau);
    # Score statistic
    Ts = qForm(K=K,s=S);
    # Calculate p-value
    p = pchisq(q=Ts,df=1,lower.tail=F);
    # Results
    Res = c(Ts,1,p);
    names(Res) = c("Score","df","p");
    return(Res);
  }
  return(Out)
}