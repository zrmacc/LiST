// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// fitNorm
SEXP fitNorm(const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::MatrixXd> Z);
RcppExport SEXP _NST_fitNorm(SEXP ySEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(fitNorm(y, Z));
    return rcpp_result_gen;
END_RCPP
}
// tr
SEXP tr(const Eigen::Map<Eigen::MatrixXd> A);
RcppExport SEXP _NST_tr(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(tr(A));
    return rcpp_result_gen;
END_RCPP
}
// fastMMp
SEXP fastMMp(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _NST_fastMMp(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(fastMMp(A, B));
    return rcpp_result_gen;
END_RCPP
}
// fastT
SEXP fastT(const Eigen::Map<Eigen::MatrixXd> A);
RcppExport SEXP _NST_fastT(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(fastT(A));
    return rcpp_result_gen;
END_RCPP
}
// fastIP
SEXP fastIP(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _NST_fastIP(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(fastIP(A, B));
    return rcpp_result_gen;
END_RCPP
}
// fastInv
SEXP fastInv(const Eigen::Map<Eigen::MatrixXd> A);
RcppExport SEXP _NST_fastInv(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(fastInv(A));
    return rcpp_result_gen;
END_RCPP
}
// fastDet
SEXP fastDet(const Eigen::Map<Eigen::MatrixXd> A);
RcppExport SEXP _NST_fastDet(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(fastDet(A));
    return rcpp_result_gen;
END_RCPP
}
// fastQF
SEXP fastQF(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> A);
RcppExport SEXP _NST_fastQF(SEXP XSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(fastQF(X, A));
    return rcpp_result_gen;
END_RCPP
}
// SchurC
SEXP SchurC(const Eigen::Map<Eigen::MatrixXd> I11, const Eigen::Map<Eigen::MatrixXd> I22, const Eigen::Map<Eigen::MatrixXd> I12);
RcppExport SEXP _NST_SchurC(SEXP I11SEXP, SEXP I22SEXP, SEXP I12SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type I11(I11SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type I22(I22SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type I12(I12SEXP);
    rcpp_result_gen = Rcpp::wrap(SchurC(I11, I22, I12));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NST_fitNorm", (DL_FUNC) &_NST_fitNorm, 2},
    {"_NST_tr", (DL_FUNC) &_NST_tr, 1},
    {"_NST_fastMMp", (DL_FUNC) &_NST_fastMMp, 2},
    {"_NST_fastT", (DL_FUNC) &_NST_fastT, 1},
    {"_NST_fastIP", (DL_FUNC) &_NST_fastIP, 2},
    {"_NST_fastInv", (DL_FUNC) &_NST_fastInv, 1},
    {"_NST_fastDet", (DL_FUNC) &_NST_fastDet, 1},
    {"_NST_fastQF", (DL_FUNC) &_NST_fastQF, 2},
    {"_NST_SchurC", (DL_FUNC) &_NST_SchurC, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_NST(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
