# Purpose: Score test for normal linear models
# Updated: 180719

#' @useDynLib NST
#' @importFrom Rcpp sourceCpp
NULL

#' Score Test for Normal Linear Models
#' 
#' Performs a score test of \eqn{H_{0}:\beta_{1}=\beta_{10}}. The score test is 
#' specified using a logical vector \code{L}, with as many entries as columns in
#' the model matrix \code{Z}. The values of \code{L} set to \code{T} are fixed
#' at \eqn{\beta_{10}} under the null. The values of \code{L} set to \code{F}
#' are estimated.
#' 
#' 
#' @param y Numeric response vector.
#' @param Z Numeric model matrix.
#' @param L Logical vector, with as many entires as columns in the model matrix,
#'   indicating which columns have fixed coefficients under the null.
#' @param b10 Value of the regression coefficient for the selected columns under 
#'   the null. Defaults to zero.
#' @return A numeric vector containing the score statistic, the degrees of
#'   freedom, and a p-value estimated using the chi-square distribution.
#'   
#' @importFrom stats model.matrix pchisq
#' @export
#' 
#' @examples
#' y = D$y;
#' Z = D$Z;
#' # Overall test for covariate effects, should reject.
#' Score.nlm(y=y,Z=Z,L=c(FALSE,TRUE,TRUE,TRUE));
#' # Individual test of b3=0, should not reject.
#' Score.nlm(y=y,Z=Z,L=c(FALSE,FALSE,FALSE,TRUE));
#' # Joint test for the effects of b1=0 and b2=0, should reject.
#' Score.nlm(y=y,Z=Z,L=c(FALSE,TRUE,TRUE,FALSE));
#' # Individual test of b1=0.1, should not reject.
#' Score.nlm(y=y,Z=Z,L=c(FALSE,TRUE,FALSE,FALSE),b10=c(0.1));

Score.nlm = function(y,Z,L,b10){
  ## Input checks
  if(!is.vector(y)){stop("A numeric vector is expected for y.")};
  if(!is.matrix(Z)){stop("A numeric matrix is expected for Z.")};
  if(!is.logical(L)){stop("A logical vector is expected for L.")};
  if(length(L)!=ncol(Z)){stop("L should have as many entries as columns in Z.")};
  # Degrees of freedom
  df = sum(L);
  if(df==0){stop("At least 1 entry of L should be TRUE.")}
  if(df==length(L)){stop("At least 1 entry of L should be FALSE.")};
  # Check for missingness
  Miss = sum(is.na(y))+sum(is.na(Z));
  if(Miss>0){stop("None of ys, Zt, or Zs, should contain missing values.")};
  # Null coefficient
  if(missing(b10)){b10=rep(0,times=df)};
  
  # Partition target design
  # Zt1 is fixed under the null.
  # Zt2 is estimated under the null.
  Z1 = Z[,L,drop=F];
  Z2 = Z[,!L,drop=F];
  # Adjust response for fixed component
  y = y-as.numeric(fastMMp(Z1,b10));
  # Fit the null model
  M0 = fitNorm(y=y,Z=Z2);
  # Calculate the score
  U = fastIP(Z1,M0$Resid);
  # Effincient information
  I11 = fastIP(Z1,Z1);
  I12 = fastIP(Z1,Z2);
  I22 = (M0$Ibb)*M0$Tau;
  # Note: I22 = fastIP(Z2,Z2);
  V = fastInv(SchurC(I11=I11,I22=I22,I12=I12));
  # Score statistic
  Ts = as.numeric(fastQF(X=U,A=V)/M0$Tau);
  # Calculate p-value
  p = pchisq(q=Ts,df=df,lower.tail=F);
  # Output
  Out = c(Ts,df,p);
  names(Out) = c("Score","df","p");
  return(Out);
}