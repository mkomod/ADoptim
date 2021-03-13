#ifndef AD_OPTIM_H
#define AD_OPTIM_H

#include <string>
#include <Rcpp.h>
#include <R/R_ext//Applic.h>

#include "adfunc.hpp"


namespace adoptim {

template <typename T>
Rcpp::List optim(T &func, std::vector<double> x, 
    std::string method = "Nelder-Mead", int type = 1, int maxit = 10000, 
    double fnscale = 1.0, int trace = 0, double abstol = R_NegInf, 
    double reltol = 1.5e-8, double alpha = 1.0, double beta = 0.5, 
    double gamma = 2.0, int report = 10, double *parscale = nullptr)
{
    double Fmin;
    int fail;
    int fncount = 0;
    int grcount = 0;
    int n = x.size();
    double x_init[n];
    double x_out[x.size()];
    double pscale[n];

    if (method != "Nelder-Mead" && method != "BFGS" && method != "CG")
	Rcpp::stop("adoptim::optim(): unknown 'method'");

    for (int i = 0; i < n; ++i) {
	x_init[i] = x.at(i);
	x_out[i] = 0;
    }

    if (parscale == nullptr) {
	for (int i = 0; i < n; ++i)
	    pscale[i] = 1;
    } else {
	for (int i = 0; i < n; ++i)
	    pscale[i] = pscale[i];
    }

    func.os.fnscale_ = fnscale;
    func.os.parscale_.assign(pscale, pscale + n);

    if (method == "Nelder-Mead") {
	nmmin(n, x_init, x_out, &Fmin, &fminfn, &fail, abstol, reltol, &func, 
		alpha, beta, gamma, trace, &fncount, maxit);
    } else if (method == "BFGS") {
	double b[n];
	int mask[n];
	for (int i = 0; i < n; ++i) {
	    b[i] = x_init[i] / pscale[i];
	    mask[i] = 1;
	}
	vmmin(n, b, &Fmin, &fminfn, &fmingr, maxit, trace, mask, abstol, 
		reltol, report, &func, &fncount, &grcount, &fail);
    } else if (method == "CG") {
	cgmin(n, x_init, x_out, &Fmin, &fminfn, &fmingr, &fail, abstol, reltol,
		&func, type, trace, &fncount, &grcount, maxit);
    }
    
    std::vector<double> par(n, 0);
    par.assign(x_out, x_out + n);

    return Rcpp::List::create(
	    Rcpp::Named("fmin") = Fmin,
	    Rcpp::Named("par") = par,
	    Rcpp::Named("fncount") = fncount,
	    Rcpp::Named("grcount") = grcount,
	    Rcpp::Named("fail") = fail
    );
};


} // namespace adoptim

#endif
