# Purpose: Prefit Score Test
# Updated: 180920

#' Prefit Model Score Test
#' 
#' Performs a joint score test on the covariates in G using a previously fit null model \code{M0}.
#' 
#' @param M0 Prefit null model. 
#' @param G Covariates tested jointly for association with the outcome. 
#' @param X Model matrix used to fit \code{M0}.
#' @param method Either "asymptotic" or "perturbation".
#' @param B Score perturbations.
#' @param parallel Run in parallel? Must reigster parallel backend first. 
#' 
#' @importFrom foreach foreach '%dopar%' registerDoSEQ
#' @importFrom stats rbinom
#' @export
#' 
#' @examples
#' \dontrun{
#' Fit null model
#' M0 = Fit.LinReg(y=y,X=X);
#' # Matrix to test for association
#' G = matrix(rnorm(2e3),nrow=1e3);
#' # Prefit model score test
#' pScore(M0=M,G=G,X=X);
#' # Standard score test
#' Score(y=y,X=cbind(X,G),L=c(rep(F,ncol(X)),rep(T,ncol(G))));
#' }

pScore = function(M0,G,X,method="asymptotic",B=1e3,parallel=F){
  ## Input checks
  if(class(M0)!="fit"){stop("A null model of class fit is expected for M0.")}
  if(!is.matrix(G)){stop("A numeric matrix is expected for G.")};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  # Method
  Choices = c("asymptotic","perturbation");
  if(!(method%in%Choices)){stop("Select method from among asymptotic or perturbation.")};
  # Unparallelize
  if(!parallel){registerDoSEQ()};
  
  # Degrees of freedom
  k = ncol(G);
  n = nrow(G);
  # Calculate score
  e = M0@Residuals;
  U = matIP(G,e);
  # Residual variance
  V = M0@V;
  # Effincient information
  I11 = matIP(G,G);
  I12 = matIP(G,X);
  I22 = (M0@Information)*V;
  E = SchurC(I11,I22,I12);
  # Score statistic
  Ts = (1/V)*as.numeric(matQF(U,matInv(E)));
  # Estimate p-value
  if(method=="asymptotic"){
    # Asymptotic
    p = pchisq(q=Ts,df=k,lower.tail=F);
  } else {
    # Perturbations
    R = foreach(i=1:(B-1),.combine=c) %dopar% {
      # Weights
      u = 2*rbinom(n=n,size=1,prob=0.5)-1;
      # Perturbed score
      Ub = matIP(G,u*e);
      # Statistic
      Tb = (1/V)*as.numeric(matQF(Ub,matInv(E)));
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