// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

//' Generate Exchangable Correlation Structure
//' 
//' @param n Dimension
//' @param r Correlation
//' @export
// [[Rcpp::export]]
SEXP exchCorr(const int n, const double r){
  const Eigen::MatrixXd J = Eigen::MatrixXd::Constant(n, n, r);
  const Eigen::VectorXd i = Eigen::VectorXd::Constant(n,(1-r));
  const Eigen::MatrixXd I = i.asDiagonal();
  const Eigen::MatrixXd Out = I + J;
  return Rcpp::wrap(Out);
}

//' Estimate Alpha under H0
//' 
//' Estimates coefficient on nuisance parameter under the \eqn{H_{0}:\beta=0}.
//' @param X Design matrix for alpha, numeric.
//' @param R Correlation structure, numeric.
//' @param y Response, numeric.
//' 
// [[Rcpp::export]]
SEXP Alpha0(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> R, const Eigen::Map<Eigen::VectorXd> y){
  const Eigen::MatrixXd Ri = R.inverse();
  const Eigen::MatrixXd XtRi = (X.transpose()*Ri);
  const Eigen::MatrixXd XtRiXi = (XtRi*X).completeOrthogonalDecomposition().pseudoInverse();
  const Eigen::VectorXd out = XtRiXi*(XtRi*y);
  return Rcpp::wrap(out);
}

//' Estimates Tau under H1
//' 
//' Estimate the variance component \eqn{\tau} under the full model.
//' 
//' @param D Overall design matrix, numeric.
//' @param R Correlation structure, numeric.
//' @param y Response, numeric.
//' 
// [[Rcpp::export]]

SEXP Tau1(const Eigen::Map<Eigen::MatrixXd> D, const Eigen::Map<Eigen::MatrixXd> R, const Eigen::Map<Eigen::VectorXd> y){
  const Eigen::MatrixXd Ri = R.inverse();
  const Eigen::MatrixXd DtRi = D.transpose()*Ri;
  const Eigen::MatrixXd DtRiDi = (DtRi*D).completeOrthogonalDecomposition().pseudoInverse();
  const Eigen::MatrixXd yhat = D*(DtRiDi*(DtRi*y));
  const Eigen::VectorXd resid = y-yhat;
  const double SS = (resid.transpose())*Ri*resid;
  const int n = D.rows();
  const int k = D.cols();
  const double tau = SS/(n-k);
  return Rcpp::wrap(tau);
}

//' Estimate Information Matrix
//' 
//' Estimate information matrix for \eqn{(\beta,\alpha)}.
//' 
//' @param X Design matrix for alpha, numeric.
//' @param G Design matrix for beta, numeric. 
//' @param R Correlation structure, numeric.
//' @param tau Variance component, numeric.
//' 
// [[Rcpp::export]]

SEXP Info(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> G, 
          const Eigen::Map<Eigen::MatrixXd> R, const double tau){
  const Eigen::MatrixXd Ri = R.inverse();
  // Iaa
  const Eigen::MatrixXd Iaa = (1/tau)*(X.transpose()*Ri)*X;
  // Iba
  const Eigen::MatrixXd Iba = (1/tau)*(G.transpose()*Ri)*X;
  // Ibb
  const Eigen::MatrixXd Ibb = (1/tau)*(G.transpose()*Ri)*G;
  return Rcpp::List::create(Rcpp::Named("Ibb") = Ibb,Rcpp::Named("Iba") = Iba,Rcpp::Named("Iaa") = Iaa);
}

//' Schur Complement
//' 
//' Calcualtes the matrix \deqn{I_{gg}-I_{gh}I_{hh}^{-1}I_{gh}'}
//' @param Igg First information matrix.
//' @param Ihh Second information matrix.
//' @param Igh Cross information matrix.
//' @param inv Return the inverse of 
//' 
//' @export  
// [[Rcpp::export]]

SEXP SchurC(const Eigen::Map<Eigen::MatrixXd> Igg, const Eigen::Map<Eigen::MatrixXd> Ihh, 
            const Eigen::Map<Eigen::MatrixXd> Igh, const bool inv){
  // Inverse of B
  const Eigen::MatrixXd Ihhi = Ihh.completeOrthogonalDecomposition().pseudoInverse();
  // Schur Complement
  const Eigen::MatrixXd S = Igg - Igh * Ihhi * Igh.transpose();
  // Inverse Schur Complement
  if(inv==TRUE){
    const Eigen::MatrixXd Out = S.completeOrthogonalDecomposition().pseudoInverse();
    return Rcpp::wrap(Out);
  }
  return Rcpp::wrap(S);
}

//' Calculate Score under H0
//' 
//' Calculate the score vector for \eqn{\beta} under \eqn{H_{0}:\beta=0}.
//' 
//' @param e0 Residuals under H0, numeric.
//' @param G Design matrix for beta, numeric. 
//' @param R Correlation structure, numeric.
//' @param tau Variance component, numeric.
//' 
// [[Rcpp::export]]

SEXP ScoreB(const Eigen::Map<Eigen::VectorXd> e0, const Eigen::Map<Eigen::MatrixXd> G, 
          const Eigen::Map<Eigen::MatrixXd> R, const double tau){
  const Eigen::MatrixXd Ri = R.inverse();
  const Eigen::VectorXd S = (1/tau) * G.transpose() * Ri * e0;
  return Rcpp::wrap(S);
}

//' Quadratic Form
//' 
//' @param K Numeric matrix.
//' @param s Numeric vector
//' @export
// [[Rcpp::export]]
SEXP qForm(const Eigen::Map<Eigen::MatrixXd> K, const Eigen::VectorXd s){
  const double Out = (s.transpose() * K * s);
  return Rcpp::wrap(Out);
}