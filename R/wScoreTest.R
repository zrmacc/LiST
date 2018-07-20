# Purpose: Variance Component Score Test
# Updated: 180719

#' Variance Component Score Test
#' 
#' Performs a weighted score test of \eqn{H_{0}:\beta_{1}=0}. The score
#' test is specified using a logical vector \code{L}, with as many entries as
#' columns in the model matrix \code{Z}. The values of \code{L} set to \code{T}
#' are fixed at \eqn{\beta_{10}} under the null. The values of \code{L} set to
#' \code{F} are estimated. The weight vector \code{w} should contain one entry for each
#' element of \code{L} set to \code{T}.
#' 
#' @param y Response vector.
#' @param Z Numeric model matrix.
#' @param w Numeric weight vector. 
#' @param L Logical vector, with as many entires as columns in the model matrix,
#'   indicating which columns have fixed coefficients under the null.
#' @return Returns a list with the following components:
#' \item{Score}{The weighted score statistic.}
#' \item{Lambda}{The eigenvalues of \eqn{\Xi}.}
#' \item{p}{The estimated p-value.}
#' 
#' @importFrom CompQuadForm davies
#' @export
#' 
#' @examples
#' y = D$y;
#' Z = D$Z;
#' # Overall test for covariate effects, should reject.
#' wScore.nlm(y=y,Z=Z,L=c(FALSE,TRUE,TRUE,TRUE))$p;
#' # Overall test, directed against b3
#' wScore.nlm(y=y,Z=Z,L=c(FALSE,TRUE,TRUE,TRUE),w=c(1/5,1/5,3/5))$p;
#' # Joint test of b1=0 and b2=0, directed against b1
#' wScore.nlm(y=y,Z=Z,L=c(FALSE,TRUE,TRUE,FALSE),w=c(3/4,1/4))$p;
#' # Joint test of b1=0 and b3=0, directed against b1
#' wScore.nlm(y=y,Z=Z,L=c(FALSE,TRUE,FALSE,TRUE),w=c(3/4,1/4))$p;

wScore.nlm = function(y,Z,L,w){
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
  if(Miss>0){stop("Neither y nor Z should contain missing values.")};
  # Weights
  if(missing(w)){w=rep(1,times=df)};
  if(length(w)!=df){stop("A numeric vector of length sum(L) is expected for w.")};
  
  # Partition target design
  # Zt1 is fixed under the null.
  # Zt2 is estimated under the null.
  Z1 = Z[,L,drop=F];
  Z2 = Z[,!L,drop=F];
  
  # Fit the null model
  M0 = fitNorm(y=y,Z=Z2);
  # Calculate the score
  U = fastIP(Z1,M0$Resid);
  # Weight matrix
  W = diag(w);
  W2 = diag(w^2);
  # Calculate the score statistic
  Ts = as.numeric(fastQF(X=U,A=W2));
  ## Construct Xi
  I11 = fastIP(Z1,Z1);
  I12 = fastIP(Z1,Z2);
  I22 = (M0$Ibb)*M0$Tau;
  # Note: I22 = fastIP(Z2,Z2);
  V = SchurC(I11=I11,I22=I22,I12=I12);
  Xi = M0$Tau*fastQF(X=W,A=V);
  # Eigenvalues
  lambda = eigen(x=Xi,symmetric=T,only.values=T)$values;
  lambda = round(lambda,digits=6);
  # Calculate p-value
  p = CompQuadForm::davies(q=Ts,lambda=lambda)$Qq;
  # Output
  Out = list("Score"=Ts,"Eigenvalues"=lambda,"p"=p);
  return(Out);
}