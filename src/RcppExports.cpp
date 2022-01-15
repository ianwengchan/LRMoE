// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ColMaxIdx
SEXP ColMaxIdx(SEXP w);
RcppExport SEXP _LRMoE_ColMaxIdx(SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(ColMaxIdx(w));
    return rcpp_result_gen;
END_RCPP
}
// EMalphadQ
SEXP EMalphadQ(SEXP x, SEXP zj, SEXP z, SEXP p);
RcppExport SEXP _LRMoE_EMalphadQ(SEXP xSEXP, SEXP zjSEXP, SEXP zSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type zj(zjSEXP);
    Rcpp::traits::input_parameter< SEXP >::type z(zSEXP);
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(EMalphadQ(x, zj, z, p));
    return rcpp_result_gen;
END_RCPP
}
// EMalphadQ2
SEXP EMalphadQ2(SEXP x, SEXP z, SEXP p, SEXP q);
RcppExport SEXP _LRMoE_EMalphadQ2(SEXP xSEXP, SEXP zSEXP, SEXP pSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type z(zSEXP);
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< SEXP >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(EMalphadQ2(x, z, p, q));
    return rcpp_result_gen;
END_RCPP
}
// EMwwbetadQ2
SEXP EMwwbetadQ2(SEXP zj, SEXP pj, SEXP p, double betajl, SEXP wbetal, SEXP wl, SEXP tl);
RcppExport SEXP _LRMoE_EMwwbetadQ2(SEXP zjSEXP, SEXP pjSEXP, SEXP pSEXP, SEXP betajlSEXP, SEXP wbetalSEXP, SEXP wlSEXP, SEXP tlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type zj(zjSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pj(pjSEXP);
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type betajl(betajlSEXP);
    Rcpp::traits::input_parameter< SEXP >::type wbetal(wbetalSEXP);
    Rcpp::traits::input_parameter< SEXP >::type wl(wlSEXP);
    Rcpp::traits::input_parameter< SEXP >::type tl(tlSEXP);
    rcpp_result_gen = Rcpp::wrap(EMwwbetadQ2(zj, pj, p, betajl, wbetal, wl, tl));
    return rcpp_result_gen;
END_RCPP
}
// EMwwdQ
SEXP EMwwdQ(SEXP z, SEXP zmarg, SEXP p, SEXP betal, SEXP tl, SEXP wwl, float sigma);
RcppExport SEXP _LRMoE_EMwwdQ(SEXP zSEXP, SEXP zmargSEXP, SEXP pSEXP, SEXP betalSEXP, SEXP tlSEXP, SEXP wwlSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type z(zSEXP);
    Rcpp::traits::input_parameter< SEXP >::type zmarg(zmargSEXP);
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< SEXP >::type betal(betalSEXP);
    Rcpp::traits::input_parameter< SEXP >::type tl(tlSEXP);
    Rcpp::traits::input_parameter< SEXP >::type wwl(wwlSEXP);
    Rcpp::traits::input_parameter< float >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(EMwwdQ(z, zmarg, p, betal, tl, wwl, sigma));
    return rcpp_result_gen;
END_RCPP
}
// EMwwdQ2
SEXP EMwwdQ2(SEXP zmarg, SEXP p, SEXP betal, SEXP tl, double sigma);
RcppExport SEXP _LRMoE_EMwwdQ2(SEXP zmargSEXP, SEXP pSEXP, SEXP betalSEXP, SEXP tlSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type zmarg(zmargSEXP);
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    Rcpp::traits::input_parameter< SEXP >::type betal(betalSEXP);
    Rcpp::traits::input_parameter< SEXP >::type tl(tlSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(EMwwdQ2(zmarg, p, betal, tl, sigma));
    return rcpp_result_gen;
END_RCPP
}
// sumBinomialY
double sumBinomialY(double n, double p, double lower_, double upper_);
RcppExport SEXP _LRMoE_sumBinomialY(SEXP nSEXP, SEXP pSEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< double >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumBinomialY(n, p, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// sumBinomialYObs
SEXP sumBinomialYObs(SEXP n_, SEXP p_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_sumBinomialYObs(SEXP n_SEXP, SEXP p_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type n_(n_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p_(p_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumBinomialYObs(n_, p_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// sumBinomialYLat
SEXP sumBinomialYLat(SEXP n_, SEXP p_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_sumBinomialYLat(SEXP n_SEXP, SEXP p_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type n_(n_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p_(p_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumBinomialYLat(n_, p_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intBurrLogYObs
SEXP intBurrLogYObs(SEXP k_, SEXP c_, SEXP lambda_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intBurrLogYObs(SEXP k_SEXP, SEXP c_SEXP, SEXP lambda_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type k_(k_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type c_(c_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intBurrLogYObs(k_, c_, lambda_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intBurrLogYLat
SEXP intBurrLogYLat(SEXP k_, SEXP c_, SEXP lambda_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intBurrLogYLat(SEXP k_SEXP, SEXP c_SEXP, SEXP lambda_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type k_(k_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type c_(c_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intBurrLogYLat(k_, c_, lambda_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intBurrPolYObs
SEXP intBurrPolYObs(SEXP k_, SEXP c_, SEXP lambda_, SEXP cc_, SEXP ll_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intBurrPolYObs(SEXP k_SEXP, SEXP c_SEXP, SEXP lambda_SEXP, SEXP cc_SEXP, SEXP ll_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type k_(k_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type c_(c_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type cc_(cc_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type ll_(ll_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intBurrPolYObs(k_, c_, lambda_, cc_, ll_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intBurrPolYLat
SEXP intBurrPolYLat(SEXP k_, SEXP c_, SEXP lambda_, SEXP cc_, SEXP ll_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intBurrPolYLat(SEXP k_SEXP, SEXP c_SEXP, SEXP lambda_SEXP, SEXP cc_SEXP, SEXP ll_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type k_(k_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type c_(c_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type cc_(cc_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type ll_(ll_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intBurrPolYLat(k_, c_, lambda_, cc_, ll_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intGammaLogYObs
SEXP intGammaLogYObs(SEXP m_, SEXP theta_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intGammaLogYObs(SEXP m_SEXP, SEXP theta_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type m_(m_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intGammaLogYObs(m_, theta_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intGammaLogYLat
SEXP intGammaLogYLat(SEXP m_, SEXP theta_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intGammaLogYLat(SEXP m_SEXP, SEXP theta_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type m_(m_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intGammaLogYLat(m_, theta_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intInvGaussLogYObs
SEXP intInvGaussLogYObs(SEXP mu_, SEXP lambda_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intInvGaussLogYObs(SEXP mu_SEXP, SEXP lambda_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intInvGaussLogYObs(mu_, lambda_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intInvGaussLogYLat
SEXP intInvGaussLogYLat(SEXP mu_, SEXP lambda_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intInvGaussLogYLat(SEXP mu_SEXP, SEXP lambda_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intInvGaussLogYLat(mu_, lambda_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intInvGaussYObs
SEXP intInvGaussYObs(SEXP mu_, SEXP lambda_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intInvGaussYObs(SEXP mu_SEXP, SEXP lambda_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intInvGaussYObs(mu_, lambda_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intInvGaussYLat
SEXP intInvGaussYLat(SEXP mu_, SEXP lambda_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intInvGaussYLat(SEXP mu_SEXP, SEXP lambda_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intInvGaussYLat(mu_, lambda_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intInvGaussInvYObs
SEXP intInvGaussInvYObs(SEXP mu_, SEXP lambda_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intInvGaussInvYObs(SEXP mu_SEXP, SEXP lambda_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intInvGaussInvYObs(mu_, lambda_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intInvGaussInvYLat
SEXP intInvGaussInvYLat(SEXP mu_, SEXP lambda_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intInvGaussInvYLat(SEXP mu_SEXP, SEXP lambda_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intInvGaussInvYLat(mu_, lambda_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// sumNegativeBinomialY
double sumNegativeBinomialY(double n, double p, double lower_, double upper_);
RcppExport SEXP _LRMoE_sumNegativeBinomialY(SEXP nSEXP, SEXP pSEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< double >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumNegativeBinomialY(n, p, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// sumNegativeBinomialYObs
SEXP sumNegativeBinomialYObs(SEXP n_, SEXP p_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_sumNegativeBinomialYObs(SEXP n_SEXP, SEXP p_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type n_(n_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p_(p_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumNegativeBinomialYObs(n_, p_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// sumNegativeBinomialYLat
SEXP sumNegativeBinomialYLat(SEXP n_, SEXP p_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_sumNegativeBinomialYLat(SEXP n_SEXP, SEXP p_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type n_(n_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p_(p_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumNegativeBinomialYLat(n_, p_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// sumNegativeBinomialLfacY
double sumNegativeBinomialLfacY(double n, double p, double nn, double lower_, double upper_);
RcppExport SEXP _LRMoE_sumNegativeBinomialLfacY(SEXP nSEXP, SEXP pSEXP, SEXP nnSEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< double >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< double >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumNegativeBinomialLfacY(n, p, nn, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// sumNegativeBinomialLfacYObs
SEXP sumNegativeBinomialLfacYObs(SEXP n_, SEXP p_, SEXP nn_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_sumNegativeBinomialLfacYObs(SEXP n_SEXP, SEXP p_SEXP, SEXP nn_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type n_(n_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p_(p_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type nn_(nn_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumNegativeBinomialLfacYObs(n_, p_, nn_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// sumNegativeBinomialLfacYLat
SEXP sumNegativeBinomialLfacYLat(SEXP n_, SEXP p_, SEXP nn_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_sumNegativeBinomialLfacYLat(SEXP n_SEXP, SEXP p_SEXP, SEXP nn_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type n_(n_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p_(p_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type nn_(nn_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumNegativeBinomialLfacYLat(n_, p_, nn_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// sumPoissonY
double sumPoissonY(double mu, double lower_, double upper_);
RcppExport SEXP _LRMoE_sumPoissonY(SEXP muSEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< double >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumPoissonY(mu, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// sumPoissonYObs
SEXP sumPoissonYObs(SEXP mu_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_sumPoissonYObs(SEXP mu_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumPoissonYObs(mu_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// sumPoissonYLat
SEXP sumPoissonYLat(SEXP mu_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_sumPoissonYLat(SEXP mu_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumPoissonYLat(mu_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intWeibullLogYObs
SEXP intWeibullLogYObs(SEXP k_, SEXP lambda_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intWeibullLogYObs(SEXP k_SEXP, SEXP lambda_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type k_(k_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intWeibullLogYObs(k_, lambda_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intWeibullLogYLat
SEXP intWeibullLogYLat(SEXP k_, SEXP lambda_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intWeibullLogYLat(SEXP k_SEXP, SEXP lambda_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type k_(k_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intWeibullLogYLat(k_, lambda_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intWeibullPowYObs
SEXP intWeibullPowYObs(SEXP k_, SEXP lambda_, SEXP p_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intWeibullPowYObs(SEXP k_SEXP, SEXP lambda_SEXP, SEXP p_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type k_(k_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p_(p_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intWeibullPowYObs(k_, lambda_, p_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// intWeibullPowYLat
SEXP intWeibullPowYLat(SEXP k_, SEXP lambda_, SEXP p_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_intWeibullPowYLat(SEXP k_SEXP, SEXP lambda_SEXP, SEXP p_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type k_(k_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p_(p_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(intWeibullPowYLat(k_, lambda_, p_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// sumZTPoissonY
double sumZTPoissonY(double mu, double lower_, double upper_);
RcppExport SEXP _LRMoE_sumZTPoissonY(SEXP muSEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< double >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumZTPoissonY(mu, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// sumZTPoissonYObs
SEXP sumZTPoissonYObs(SEXP mu_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_sumZTPoissonYObs(SEXP mu_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumZTPoissonYObs(mu_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// sumZTPoissonYLat
SEXP sumZTPoissonYLat(SEXP mu_, SEXP lower_, SEXP upper_);
RcppExport SEXP _LRMoE_sumZTPoissonYLat(SEXP mu_SEXP, SEXP lower_SEXP, SEXP upper_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_(lower_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_(upper_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumZTPoissonYLat(mu_, lower_, upper_));
    return rcpp_result_gen;
END_RCPP
}
// XAPlusYZB
SEXP XAPlusYZB(SEXP x, SEXP a, SEXP y, SEXP z, SEXP b);
RcppExport SEXP _LRMoE_XAPlusYZB(SEXP xSEXP, SEXP aSEXP, SEXP ySEXP, SEXP zSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type a(aSEXP);
    Rcpp::traits::input_parameter< SEXP >::type y(ySEXP);
    Rcpp::traits::input_parameter< SEXP >::type z(zSEXP);
    Rcpp::traits::input_parameter< SEXP >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(XAPlusYZB(x, a, y, z, b));
    return rcpp_result_gen;
END_RCPP
}
// XColMinusY
SEXP XColMinusY(SEXP x, SEXP y);
RcppExport SEXP _LRMoE_XColMinusY(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(XColMinusY(x, y));
    return rcpp_result_gen;
END_RCPP
}
// XPlusYColTimesZ
SEXP XPlusYColTimesZ(SEXP x, SEXP y, SEXP z);
RcppExport SEXP _LRMoE_XPlusYColTimesZ(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type y(ySEXP);
    Rcpp::traits::input_parameter< SEXP >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(XPlusYColTimesZ(x, y, z));
    return rcpp_result_gen;
END_RCPP
}
// XPlusYZ
SEXP XPlusYZ(SEXP x, SEXP y, SEXP z);
RcppExport SEXP _LRMoE_XPlusYZ(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type y(ySEXP);
    Rcpp::traits::input_parameter< SEXP >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(XPlusYZ(x, y, z));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LRMoE_ColMaxIdx", (DL_FUNC) &_LRMoE_ColMaxIdx, 1},
    {"_LRMoE_EMalphadQ", (DL_FUNC) &_LRMoE_EMalphadQ, 4},
    {"_LRMoE_EMalphadQ2", (DL_FUNC) &_LRMoE_EMalphadQ2, 4},
    {"_LRMoE_EMwwbetadQ2", (DL_FUNC) &_LRMoE_EMwwbetadQ2, 7},
    {"_LRMoE_EMwwdQ", (DL_FUNC) &_LRMoE_EMwwdQ, 7},
    {"_LRMoE_EMwwdQ2", (DL_FUNC) &_LRMoE_EMwwdQ2, 5},
    {"_LRMoE_sumBinomialY", (DL_FUNC) &_LRMoE_sumBinomialY, 4},
    {"_LRMoE_sumBinomialYObs", (DL_FUNC) &_LRMoE_sumBinomialYObs, 4},
    {"_LRMoE_sumBinomialYLat", (DL_FUNC) &_LRMoE_sumBinomialYLat, 4},
    {"_LRMoE_intBurrLogYObs", (DL_FUNC) &_LRMoE_intBurrLogYObs, 5},
    {"_LRMoE_intBurrLogYLat", (DL_FUNC) &_LRMoE_intBurrLogYLat, 5},
    {"_LRMoE_intBurrPolYObs", (DL_FUNC) &_LRMoE_intBurrPolYObs, 7},
    {"_LRMoE_intBurrPolYLat", (DL_FUNC) &_LRMoE_intBurrPolYLat, 7},
    {"_LRMoE_intGammaLogYObs", (DL_FUNC) &_LRMoE_intGammaLogYObs, 4},
    {"_LRMoE_intGammaLogYLat", (DL_FUNC) &_LRMoE_intGammaLogYLat, 4},
    {"_LRMoE_intInvGaussLogYObs", (DL_FUNC) &_LRMoE_intInvGaussLogYObs, 4},
    {"_LRMoE_intInvGaussLogYLat", (DL_FUNC) &_LRMoE_intInvGaussLogYLat, 4},
    {"_LRMoE_intInvGaussYObs", (DL_FUNC) &_LRMoE_intInvGaussYObs, 4},
    {"_LRMoE_intInvGaussYLat", (DL_FUNC) &_LRMoE_intInvGaussYLat, 4},
    {"_LRMoE_intInvGaussInvYObs", (DL_FUNC) &_LRMoE_intInvGaussInvYObs, 4},
    {"_LRMoE_intInvGaussInvYLat", (DL_FUNC) &_LRMoE_intInvGaussInvYLat, 4},
    {"_LRMoE_sumNegativeBinomialY", (DL_FUNC) &_LRMoE_sumNegativeBinomialY, 4},
    {"_LRMoE_sumNegativeBinomialYObs", (DL_FUNC) &_LRMoE_sumNegativeBinomialYObs, 4},
    {"_LRMoE_sumNegativeBinomialYLat", (DL_FUNC) &_LRMoE_sumNegativeBinomialYLat, 4},
    {"_LRMoE_sumNegativeBinomialLfacY", (DL_FUNC) &_LRMoE_sumNegativeBinomialLfacY, 5},
    {"_LRMoE_sumNegativeBinomialLfacYObs", (DL_FUNC) &_LRMoE_sumNegativeBinomialLfacYObs, 5},
    {"_LRMoE_sumNegativeBinomialLfacYLat", (DL_FUNC) &_LRMoE_sumNegativeBinomialLfacYLat, 5},
    {"_LRMoE_sumPoissonY", (DL_FUNC) &_LRMoE_sumPoissonY, 3},
    {"_LRMoE_sumPoissonYObs", (DL_FUNC) &_LRMoE_sumPoissonYObs, 3},
    {"_LRMoE_sumPoissonYLat", (DL_FUNC) &_LRMoE_sumPoissonYLat, 3},
    {"_LRMoE_intWeibullLogYObs", (DL_FUNC) &_LRMoE_intWeibullLogYObs, 4},
    {"_LRMoE_intWeibullLogYLat", (DL_FUNC) &_LRMoE_intWeibullLogYLat, 4},
    {"_LRMoE_intWeibullPowYObs", (DL_FUNC) &_LRMoE_intWeibullPowYObs, 5},
    {"_LRMoE_intWeibullPowYLat", (DL_FUNC) &_LRMoE_intWeibullPowYLat, 5},
    {"_LRMoE_sumZTPoissonY", (DL_FUNC) &_LRMoE_sumZTPoissonY, 3},
    {"_LRMoE_sumZTPoissonYObs", (DL_FUNC) &_LRMoE_sumZTPoissonYObs, 3},
    {"_LRMoE_sumZTPoissonYLat", (DL_FUNC) &_LRMoE_sumZTPoissonYLat, 3},
    {"_LRMoE_XAPlusYZB", (DL_FUNC) &_LRMoE_XAPlusYZB, 5},
    {"_LRMoE_XColMinusY", (DL_FUNC) &_LRMoE_XColMinusY, 2},
    {"_LRMoE_XPlusYColTimesZ", (DL_FUNC) &_LRMoE_XPlusYColTimesZ, 3},
    {"_LRMoE_XPlusYZ", (DL_FUNC) &_LRMoE_XPlusYZ, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_LRMoE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
