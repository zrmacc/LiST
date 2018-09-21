// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//' Ordinary Least Squares
//' 
//' Fits the standard OLS model.
//' 
//' @param y Numeric vector.
//' @param X Numeric matrix.
//' 
//' @return List containing the following:
//' \item{Beta}{Regression coefficient.}
//' \item{V}{Outcome variance.}
//' \item{Ibb}{Information matrix for beta.}
//' \item{Resid}{Outcome residuals.}
//' 
// [[Rcpp::export]]

SEXP fitOLS(const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::MatrixXd> X){
  // Observations
  const int n = y.size();
  // Estimated parameters
  const int p = X.cols();
  // Information
  const Eigen::MatrixXd A = X.transpose()*X;
  // Estimate beta
  const Eigen::VectorXd b = A.ldlt().solve(X.transpose()*y);
  // Calculate residuals
  const Eigen::VectorXd eps = (y-X*b);
  // Scale
  const double qf = (eps.transpose()*eps);
  const double v = qf/(n-p);
  // Information
  const Eigen::MatrixXd Ibb = A/v;
  return Rcpp::List::create(Rcpp::Named("Beta")=b,Rcpp::Named("V")=v,Rcpp::Named("Ibb")=Ibb,Rcpp::Named("Resid")=eps);
}

//' Weighted Least Squares
//' 
//' Fits the following weighted least squares model: 
//' \eqn{y_{i}=x_{i}'\beta+\epsilon_{i}}. Here, the subject-specific residual is
//' normally distributed with mean zero and variance \eqn{\sigma^{2}/w_{i}}.
//' \eqn{w_{i}} is a known, subject-specific weight, and \eqn{\sigma} is a
//' common scale parameter.
//' 
//' @param y Response vector.
//' @param X Design matrix.
//' @param w Weight vector.
//' 
//' @return List containing the following:
//' \item{Beta}{Regression coefficient.}
//' \item{V}{Outcome variance.}
//' \item{Ibb}{Information matrix for beta.}
//' \item{Resid}{Outcome residuals.}
//'
// [[Rcpp::export]]
SEXP fitWLS(const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::MatrixXd> X, 
            const Eigen::Map<Eigen::VectorXd> w){
  // Observations
  const int n = y.size();
  // Estimated parameters
  const int p = X.cols();
  // Calculate A=X'WX 
  Eigen::MatrixXd A = Eigen::MatrixXd::Constant(p,p,0);
  for(int i=0; i<n; i++){
    A += X.row(i).transpose()*w(i)*X.row(i);
  };
  // Calculate u=X'Wy
  Eigen::VectorXd u = Eigen::VectorXd::Constant(p,0);
  for(int i=0; i<n; i++){
    u += X.row(i).transpose()*w(i)*y(i);
  };
  // Estimate beta
  const Eigen::VectorXd b = A.ldlt().solve(u);
  // Calculate residuals
  const Eigen::VectorXd eps = (y-X*b);
  // Scale
  double qf = 0;
  for(int i=0; i<n; i++){
    qf += eps(i)*w(i)*w(i)*eps(i);
  }
  const double v = qf/(n-p);
  // Information
  const Eigen::MatrixXd Ibb = A/v;
  return Rcpp::List::create(Rcpp::Named("Beta")=b,Rcpp::Named("V")=v,Rcpp::Named("Ibb")=Ibb,Rcpp::Named("Resid")=eps);
}