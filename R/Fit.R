# Purpose: Master fitting function
# Updated: 180914

#' Fit Linear Regression Model
#' 
#' @param y Numeric outcome vector.
#' @param X Numeric model matrix. Include an intercept. 
#' @param w Optional vector of known, observation-specific weights. 
#' @param sig Significance level, for CIs.
#' 
#' @importFrom methods new
#' @importFrom stats pnorm qnorm
#' @export 
#' 
#' @examples
#' \dontrun{
#' # Generate data
#' set.seed(100);
#' # Design
#' X = cbind(1,matrix(rnorm(4e3),nrow=1e3));
#' colnames(X) = c("int",paste0("x",seq(1:4)));
#' b = c(1,-0.1,0.2,-0.1,0);
#' # Linear predictor
#' h = as.numeric(MMP(X,b));
#' # Outcome
#' y = h+rnorm(1e3);
#' # Model fit
#' M = Fit.LinReg(y=y,X=X);
#' show(M);
#' }

Fit.LinReg = function(y,X=NULL,w=NULL,sig=0.05){
  # Intput check
  n = length(y);
  if(!is.numeric(y)){stop("A numeric vector is expected for y.")};
  if(is.null(X)){X = array(1,dim=c(n,1)); colnames(X)="int";};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  
  # Model fit
  if(is.null(w)){
    M = fitOLS(y=y,X=X);
  } else {
    M = fitWLS(y=y,X=X,w=w);
  }
  
  # Coefficient
  b = M$Beta;
  # Linear predictions
  h = as.numeric(MMP(X,b));
  # Information
  J = M$Ibb;
  Ji = matInv(J);
  se = sqrt(diag(Ji));
  
  # Coefficient frame
  if(is.null(colnames(X))){colnames(X)=paste0("x",seq(1:ncol(X)))};
  B = data.frame(colnames(X),b,se);
  colnames(B) = c("Coeff","Point","SE");
  rownames(J) = colnames(J) = colnames(X);
  
  # CIs
  z = qnorm(1-sig/2);
  B$L = B$Point-z*B$SE;
  B$U = B$Point+z*B$SE;
  B$p = 2*pnorm(abs(B$Point/B$SE),lower.tail=F);
  
  # Output
  Out = new(Class="fit",Coefficients=B,Eta=h,Information=J,Residuals=M$Resid,V=M$V);
  return(Out);
}