#ifndef AD_FUNC_H
#define AD_FUNC_H

#include <vector>

#include <R/R_ext/Applic.h>
#include "cppad/cppad.hpp"


namespace adoptim {

typedef CppAD::AD<double> ad_double;
typedef CPPAD_TESTVECTOR(CppAD::AD<double>) ad_vec_double;
typedef CppAD::ADFun<double> ad_fun;


struct OptStruct {
    bool has_grad_ = false;
    bool has_hess_ = false;
    double fnscale_ = 1.0;  		// scaling for objective
    std::vector<double> parscale_;    	// scaling for parameters
    int usebounds_ = 0;
};


class AD_func {
    private:
	ad_fun gr_;
    
    public:
	OptStruct os;

	virtual ad_double fn(const ad_vec_double &x) {
	    return x[0];
	}
	
	template <class T>
	T gr(const T &x) {
	    return gr_.Jacobian(x);
	}

	void set_grad(const std::vector<double> &x) {
	    ad_vec_double ax;
	    ax.push_vector(x);
	    CppAD::Independent(ax);

	    ad_vec_double ay(1);
	    ay[0] = fn(ax);
	    gr_.Dependent(ax, ay);

	    gr_.optimize();
	}

};


inline double fminfn(int n, double *x, void *ex) 
{
    OptStruct os(static_cast<AD_func *>(ex)->os);
    std::vector<double> x_(n, 0);

    x_.assign(x, x+n);
    ad_vec_double ax_;
    ax_.push_vector(x_);
    
    for (int i = 0; i < n; ++i)
	x[i] /= os.parscale_.at(i);

    return CppAD::Value(static_cast<AD_func *>(ex)->fn(ax_)) / os.fnscale_;
}


inline void fmingr(int n, double *x, double *gr, void *ex)
{
    OptStruct os(static_cast<AD_func *>(ex)->os);
    std::vector<double> x_(n, 0);
    x_.assign(x, x+n);

    for (int i = 0; i < n; ++i)
	x[i] /= os.parscale_.at(i);

    std::vector<double> xgr = static_cast<AD_func *>(ex)->gr(x_);
    for (auto i = 0; i != n; ++i)
	gr[i] = xgr.at(i) * (os.parscale_.at(i) / os.fnscale_);
}

} // namespace adoptim

#endif
