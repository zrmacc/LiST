#' @useDynLib ScoreTest
#' @importFrom Rcpp sourceCpp
NULL

#' Linear Model Score Test
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
#' @return Vector include the score statistic, the degrees of freedom, and a
#'   p-value based on the chi-square distribution.
#' 
#' @importFrom stats model.matrix pchisq
#' @export
#' 
#' @examples 
#' ScoreTest::lmScoreTest(G=ScoreTest::G[1,],X=ScoreTest::X,y=ScoreTest::Y);

lmScoreTest = function(G,X,y,R){
  # Input formatting
  X = data.frame(X);
  G = data.frame(G);
  # Observations
  n = nrow(X);
  # Number of covariates
  k = ncol(G);
  q = ncol(X);
  # Check for correlation structure
  if(missing(R)){R=diag(n)}
  # Estimate alpha under H0
  Dx = model.matrix(~.,data=X);
  alpha0 = Alpha0(X=Dx,R=R,y=y);
  # Null residuals
  e0 = (y - Dx %*% alpha0);
  # Estimate tau under H1
  D = model.matrix(~.,data=data.frame(X,G))
  tau1 = Tau1(D=D,R=R,y=y);
  # Information matrix for (beta,alpha)
  Dg = model.matrix(~0+.,data=G);
  J = Info(X=Dx,G=Dg,R=R,tau=tau1);
  # Efficient information for beta
  K = SchurC(Igg=J$Ibb,Ihh=J$Iaa,Igh=J$Iba,inv=T);
  # Score for beta under H0
  S = ScoreB(e0=e0,G=Dg,R=R,tau=tau1);
  # Score statistic
  Ts = qForm(K=K,s=S);
  # Calculate p-value
  p = pchisq(q=Ts,df=k,lower.tail=F);
  # Output
  Out = c(Ts,k,p);
  names(Out) = c("Score","df","p");
  return(Out);
}
