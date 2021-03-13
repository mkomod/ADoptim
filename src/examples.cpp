#include <vector>
#include <Rcpp.h>

#include "adoptim.hpp"

using namespace adoptim;

class F : public AD_func {
    ad_double fn(const ad_vec_double &x) {
	return 100 * pow(x[1] - x[0]*x[0], 2) + pow(1-x[0], 2);
    }
};


Rcpp::List example()
{
    std::vector<double> x = {-1.2, 1};
    F func;
    func.set_grad(x);
    return optim(func, x);
}
