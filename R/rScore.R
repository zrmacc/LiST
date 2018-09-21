# Purpose: Repeated Score Test
# Updated: 180920

#' Repeated Score Test 
#' 
#' Individually tests each column of \code{G} for association with \code{y} 
#' while adjusting for \code{X}. \code{G} may contain missing elements, although
#' the remaining model matrices may not. \code{rScore} accelerates association
#' testing by recycling the same null model for each hypothesis test.
#' 
#' @param y Numeric response vector.
#' @param G Numeric matrix of covariates for the target outcome whose regression
#'   coefficients are zero under the null.
#' @param X Numeric model matrix included in the null model.
#' @param parallel Run association tests if parallel? Must register parallel 
#'   backend first.
#' @return Vector of asymptotic p-values, one for each column of \code{G}.
#'   
#' @importFrom plyr aaply
#' @importFrom stats pchisq
#' @export
#' 
#' @examples
#' \dontrun{
#' # Genotypes
#' G = replicate(rbinom(n=1e3,size=2,prob=0.25),n=2000);
#' G[sample(length(G),size=floor(0.01*length(G)),replace=F)] = NA;
#' storage.mode(G) = "numeric";
#' doMC::registerDoMC(cores=2);
#' # Repeated Score Test
#' P = rScore(y=y,G=G,X=X,parallel=T);
#' }

rScore = function(y,G,X,parallel=F){
  ## Input checks
  if(!is.vector(y)){stop("A numeric verctor is expected for y.")};
  if(!is.matrix(G)){stop("A numeric matrix is expected for G.")};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  # Missingness
  Miss = sum(is.na(y))+sum(is.na(X));
  if(Miss>0){stop("Neither y nor X should contain missing values.")};
  
  # Fit null model
  M0 = fitOLS(y=y,X=X);
  # Extract residuals
  e = M0$Resid;
  # Residual variance
  V = M0$V;
  # Base information matrix
  J22 = M0$Ibb*V;
  # Function to calculate score statistics
  aux = function(g){
    # Adjust for missingness
    key = !is.na(g);
    g0 = g[key];
    X0 = X[key,,drop=F];
    X1 = X[!key,,drop=F];
    e0 = e[key];
    # Information components
    I11 = matIP(g0,g0);
    I12 = matIP(g0,X0);
    I22 = J22-matIP(X1,X1);
    # Variance
    E = as.numeric(SchurC(I11,I22,I12));
    # Score
    u = as.numeric(matIP(g0,e0));
    # Test statistic
    Ts = u^2/(V*E);
    return(Ts);
  }
  # Calculate score statistics
  U = aaply(.data=G,.margins=2,.fun=aux,.parallel=parallel);
  # Calculate p-value
  p = pchisq(q=U,df=1,lower.tail=F);
  return(p);
}