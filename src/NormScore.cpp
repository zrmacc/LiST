// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

//' Normal Score Statistic
//' 
//' @param y Outcome.
//' @param X1 Model matrix whose coefficient is fixed.
//' @param b Values at which the coefficients of X1 are fixed.
//' @param X2 Model matrix whose coefficient is estimated.
//' @param estT Logical indicating tau should be estimated. If false, provide the
//'   known value.
//' @param t Variance component, if known.
//' @param useK Logical indicating that a known covariance structure is supplied.
//' @param K Covariance structure, if known.
// [[Rcpp::export]]

SEXP normScore(const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::MatrixXd> X1,
               const Eigen::Map<Eigen::VectorXd> b, const Eigen::Map<Eigen::MatrixXd> X2,
               const bool estT, const double t, const bool useK, 
               const Eigen::Map<Eigen::MatrixXd> K){
  // Observations
  const int n = y.size();
  // Estimated parameters
  const int p = X2.cols();
  // Declare Qk
  Eigen::MatrixXd Qk(n,n);
  // Construct Qk
  if(useK){
    const Eigen::MatrixXd Ki = K.completeOrthogonalDecomposition().pseudoInverse();
    Qk = Ki-X2*(X2.transpose()*K*X2).llt().solve(X2.transpose());
  } else {
    const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n,n);
    Qk = I-X2*(X2.transpose()*X2).llt().solve(X2.transpose());
  }
  // Declare tau
  double tau;
  // Construct tau
  if(estT){
    const double qf = y.transpose()*Qk*y;
    tau = qf/(n-p);
  } else {
    tau = t;
  }
  // Score
  const Eigen::VectorXd u = X1.transpose()*Qk*(y-X1*b);
  // Test statistic
  const double T0 = u.transpose()*(X1.transpose()*Qk*X1).llt().solve(u);
  const double T1 = T0/tau;
  return Rcpp::wrap(T1);
}

//' Weighted Normal Score Statistic
//' 
//' @param y Outcome.
//' @param X1 Model matrix whose coefficient is fixed.
//' @param b Values at which the coefficients of X1 are fixed.
//' @param X2 Model matrix whose coefficient is estimated.
//' @param estT Logical indicating tau should be estimated. If false, provide the
//'   known value.
//' @param t Variance component, if known.
//' @param useK Logical indicating that a known covariance structure is supplied.
//' @param K Covariance structure, if known.
//' @param W Weight matrix
// [[Rcpp::export]]

SEXP wNormScore(const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::MatrixXd> X1,
                const Eigen::Map<Eigen::VectorXd> b, const Eigen::Map<Eigen::MatrixXd> X2,
                const bool estT, const double t, const bool useK, 
                const Eigen::Map<Eigen::MatrixXd> K, const Eigen::Map<Eigen::MatrixXd> W){
  // Observations
  const int n = y.size();
  // Estimated parameters
  const int p = X2.cols();
  // Declare Qk
  Eigen::MatrixXd Qk(n,n);
  // Construct Qk
  if(useK){
    const Eigen::MatrixXd Ki = K.completeOrthogonalDecomposition().pseudoInverse();
    Qk = Ki-X2*(X2.transpose()*K*X2).llt().solve(X2.transpose());
  } else {
    const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n,n);
    Qk = I-X2*(X2.transpose()*X2).llt().solve(X2.transpose());
  }
  // Declare tau
  double tau;
  // Construct tau
  if(estT){
    const double qf = y.transpose()*Qk*y;
    tau = qf/(n-p);
  } else {
    tau = t;
  }
  // Score
  const Eigen::VectorXd u = Qk*(y-X1*b);
  // Test statistic
  const double T0 = u.transpose()*W*W*u;
  const double T1 = T0/(tau*tau);
  // Q matrix
  const Eigen::MatrixXd Q = Qk/(tau);
  return Rcpp::List::create(Rcpp::Named("Tw")=T1,Rcpp::Named("Q")=Q);
}

/////////////////////////////////////////////////

//' Normal Model
//' 
//' @param y Outcome.
//' @param Z Model matrix.
//' @param estT Logical indicating tau should be estimated. If false, provide the
//'   known value.
//' @param t Variance component, if known.
// [[Rcpp::export]]

SEXP fitNorm(const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::MatrixXd> Z,
             const bool estT, const double t){
  // Observations
  const int n = y.size();
  // Estimated parameters
  const int p = Z.cols();
  // Gram matrix
  const Eigen::MatrixXd ZtZ = Z.transpose()*Z;
  // Estimate beta
  const Eigen::VectorXd b = (ZtZ).llt().solve(Z.transpose()*y);
  // Calculate residuals
  const Eigen::VectorXd eT = (y-Z*b);
  // Declare tau
  double tau;
  if(estT){
    const double qf = (eT.transpose()*eT);
    tau = qf/(n-p);
  } else {
    tau = t;
  }
  // Information
  const Eigen::MatrixXd Ibb = ZtZ/tau;
  return Rcpp::List::create(Rcpp::Named("Beta")=b,Rcpp::Named("Tau")=tau,Rcpp::Named("Ibb")=Ibb,Rcpp::Named("eT")=eT);
}