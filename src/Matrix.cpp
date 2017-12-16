// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

// Purpose : Matrix operations

//' Matrix Inner Product
//' 
//' Forms the product \eqn{X'AX}.
//' @param X Numeric matrix.
//' @param A Numeric matrix.
// [[Rcpp::export]]
SEXP matIP(const Eigen::Map<Eigen::MatrixXd> X,const Eigen::Map<Eigen::MatrixXd> A){
  const Eigen::MatrixXd Out = (X.transpose() * A * X);
  return Rcpp::wrap(Out);
}

//' Fast Inverse
//' 
//' Forms the product \eqn{X'AX}.
//' @param A Numeric matrix.
// [[Rcpp::export]]
SEXP fastInv(const Eigen::Map<Eigen::MatrixXd> A){
  const Eigen::MatrixXd Out = A.completeOrthogonalDecomposition().pseudoInverse();
  return Rcpp::wrap(A);
}

//' Schur Complement
//' 
//' Calcualtes the matrix \deqn{I_{gg}-I_{gh}I_{hh}^{-1}I_{gh}'}
//' @param Igg First information matrix.
//' @param Ihh Second information matrix.
//' @param Igh Cross information matrix.
//' @param inv Return the inverse of the Schur complement?
//'   
// [[Rcpp::export]]
SEXP SchurC(const Eigen::Map<Eigen::MatrixXd> Igg, const Eigen::Map<Eigen::MatrixXd> Ihh, 
            const Eigen::Map<Eigen::MatrixXd> Igh, const bool inv){
  // Schur Complement
  const Eigen::MatrixXd S = Igg - Igh * Ihh.llt().solve(Igh.transpose());
  // Inverse Schur Complement
  if(inv==TRUE){
    const Eigen::MatrixXd Out = S.completeOrthogonalDecomposition().pseudoInverse();
    return Rcpp::wrap(Out);
  }
  return Rcpp::wrap(S);
}