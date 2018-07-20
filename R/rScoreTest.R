# Purpose: Repeated Score Test
# Updated: 180719

#' Repeated Score Test for Normal Models
#' 
#' Individually tests \eqn{H_{0}:\beta=0} for each column of \code{Z1}, 
#' adjusting for \code{Z2}. Missing values are permitted in \code{Z1}, 
#' though not in \code{Z2}.
#' 
#' @param y Numeric response vector.
#' @param Z1 Numeric matrix of covariates for the target outcome whose 
#'   regression coefficients are zero under the null.
#' @param Z2 Numeric matrix of covariates for the target outcome whose 
#'   regression coefficients are unconstrained, i.e. estimated, under the null.
#' @param parallel Run association tests if parallel? Must register parallel 
#'   backend first.
#' @return Vector of p-values, one for each column of \code{X1}, based on the 
#'   chi square distribution.
#'   
#' @importFrom plyr aaply
#' @importFrom stats model.matrix pchisq
#' @export
#' 
#' @examples
#' y = D$y;
#' Z = D$Z;
#' set.seed(100);
#' G = replicate(1000,rbinom(n=length(y),size=2,prob=0.25));
#' storage.mode(G) = "numeric";
#' # Test each column of G for association with y, adjusting for Z
#' p = rScore.nlm(y=y,Z1=G,Z2=Z);
#' mean(p<=0.05);

rScore.nlm = function(y,Z1,Z2,parallel=F){
  ## Input checks
  if(!is.vector(y)){stop("A numeric verctor is expected for y.")};
  if(!is.matrix(Z1)){stop("A numeric matrix is expected for Z1.")};
  if(!is.matrix(Z2)){stop("A numeric matrix is expected for Z2.")};
  # Missingness
  Miss = sum(is.na(y))+sum(is.na(Z2));
  if(Miss>0){stop("Neither y nor Z2 should contain missing values.")};
  # Fit null model
  M0 = fitNorm(y=y,Z=Z2);
  # Extract residuals
  e = M0$Resid;
  # Residual variance
  tau = M0$Tau;
  # Function to calculate score statistics
  aux = function(x){
    # Adjust for missingness
    keep = !is.na(x);
    x.obs = x[keep];
    Z.obs = Z2[keep,,drop=F];
    e.obs = e[keep];
    # Information components
    I11 = fastIP(x.obs,x.obs);
    I12 = fastIP(x.obs,Z.obs);
    I22 = fastIP(Z.obs,Z.obs);
    # Variance
    V = as.numeric(SchurC(I11=I11,I22=I22,I12=I12));
    # Score
    u = as.numeric(fastIP(x.obs,e.obs));
    # Test statistic
    Ts = u^2/(V*tau);
    return(Ts);
  }
  # Calculate score statistics
  U = aaply(.data=Z1,.margins=2,.fun=aux,.parallel=parallel);
  # Calculate p-value
  p = pchisq(q=U,df=1,lower.tail=F);
  return(p);
}