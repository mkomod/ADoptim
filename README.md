# ad-optim

Optimisation routines in R rely on functions to compute gradients. `ad-optim` is a wrapper around R's optimisations routines where function gradients are computed automatically via algorithmic differentiation (through CppAD).

## Installation

```R
R> devtools::install_github("mkomod/ad-optim")
```

## Usage 

adoptim is used via Rcpp by the inclusion of the "adoptim.hpp" header

```C++
#include "adoptim.hpp"

using namespace adoptim;

// define a function
class Func : public AD_func {
    // define the fn method
    ad_double fn(const ad_vec_double &x) {
	return 100 * pow(x[1] - x[0]*x[0], 2) + pow(1-x[0], 2);
    }
}


void some_other_function() {
    // ...
    Func func;
    std::vector<double> par = { ... }; // parmeter values
    func.set_grad(par); // initialise the gradient function

    optim(func, par); // optimise
}
```
