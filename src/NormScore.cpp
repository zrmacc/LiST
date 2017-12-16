// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

//' Estimate \eqn{\alpha} under \eqn{H_{0}}.
//' 
//' Estimates coefficient on nuisance parameter under the \eqn{H_{0}:\beta=0}.
//' @param X Design matrix for alpha, numeric.
//' @param Ri Inverse correlation structure, numeric.
//' @param y Response, numeric.
//' 
// [[Rcpp::export]]
SEXP Alpha0(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> Ri, const Eigen::Map<Eigen::VectorXd> y){
  const Eigen::MatrixXd XtRi = (X.transpose()*Ri);
  const Eigen::VectorXd a = (XtRi*X).llt().solve(XtRi*y);
  return Rcpp::wrap(a);
}

//' Estimate \eqn{tau} under \eqn{H_{0}}.
//' 
//' Estimates the variance component \eqn{\tau} under the reduced
//' model where \eqn{\beta = 0}. 
//' 
//' @param X Overall design matrix, numeric.
//' @param Ri Inverse correlation structure, numeric.
//' @param y Response, numeric.
//' @param a Estimated regression coefficient
//' 
// [[Rcpp::export]]
SEXP Tau0(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> Ri, 
          const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::VectorXd> a){
  // Residuals
  const Eigen::VectorXd e = (y-X*a);
  // Dimensions
  const int n = X.rows();
  const int k = X.cols();
  // Estimate
  const double SS = (e.transpose() * Ri * e);
  const double tau = SS/(n-k);
  return Rcpp::wrap(tau);
}

//' Score Statistic for Normal Linear Model
//' 
//' Calculates the score statistic of \eqn{H_{0}:\beta_{G}=0} in the form
//' \deqn{T_{S}=y'QG(G'QG)^{-1}G'Qy}
//' 
//' @param X Matrix whose regression coefficient is unconstrained in the null
//'   model.
//' @param G Matrix whose regression coefficient is zero in the null model.
//' @param Ri Inverse correlation structure, numeric.
//' @param y Response, numeric.
//' @param tau Estimate of tau
//'   
// [[Rcpp::export]]
SEXP nScore(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> G, 
            const Eigen::Map<Eigen::MatrixXd> Ri, const Eigen::Map<Eigen::VectorXd> y, 
            const double tau){
  // Error projection matrix
  const Eigen::MatrixXd Q = (Ri-Ri*X*(X.transpose()*Ri*X).llt().solve(X.transpose()*Ri))/tau;
  // Score statistic
  const Eigen::VectorXd s = G.transpose()*Q*y;
  const double Ts = s.transpose()*(G.transpose()*Q*G).llt().solve(s);
  return Rcpp::wrap(Ts);
}

//' Kernalized Score Statistic for Normal Linear Model
//' 
//' Calculates a kernelized score statistic of \eqn{H_{0}:\beta = 0} in the form
//' \deqn{T_{K} = \epsilon'K\epsilon} where \eqn{K = LL'}
//' 
//' @param X Matrix whose regression coefficient is unconstrained in the null
//'   model.
//' @param L Cholesky decomposition of kernel matrix.
//' @param Ri Inverse correlation structure, numeric.
//' @param y Response, numeric.
//' @param tau Estimate of tau
//'   
// [[Rcpp::export]]
SEXP knScore(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> L, 
             const Eigen::Map<Eigen::MatrixXd> Ri, const Eigen::Map<Eigen::VectorXd> y, 
             const double tau){
  // Error projection matrix
  const Eigen::MatrixXd Q = (Ri-Ri*X*(X.transpose()*Ri*X).llt().solve(X.transpose()*Ri))/tau;
  // Score statistic
  const double Tk = y.transpose()*Q*L*L.transpose()*Q*y;
  return Rcpp::List::create(Rcpp::Named("Tk")=Tk,Rcpp::Named("Q")=Q);
}