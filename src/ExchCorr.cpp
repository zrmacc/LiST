// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

// Purpose : Function to generate exchangeable correlation structure

//' Generate Exchangable Correlation Structure
//' 
//' @param n Dimension
//' @param r Correlation
//' @export
// [[Rcpp::export]]
SEXP exchCorr(const int n, const double r){
  // Generate matrix of the form r*1*t(1)
  const Eigen::MatrixXd J = Eigen::MatrixXd::Constant(n, n, r);
  // Generate diagonal matrix of the form (1-r)*I
  const Eigen::VectorXd i = Eigen::VectorXd::Constant(n,(1-r));
  const Eigen::MatrixXd I = i.asDiagonal();
  // Exchangeable correlation structure
  const Eigen::MatrixXd Out = I + J;
  return Rcpp::wrap(Out);
}