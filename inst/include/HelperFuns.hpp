using namespace R_inla;
using namespace density;
using namespace tmbutils;
using namespace Eigen;

// Manually create a sparse matrix
template<class Type>
SparseMatrix<Type> Q_spde_generalized(spde_t<Type> spde, Type kappa, int alpha=2){
  Type kappa_pow2 = kappa*kappa;
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  
  if( alpha==1 ) return kappa_pow2*spde.M0 + spde.M1;
  if( alpha==2 ) return kappa_pow4*spde.M0 + Type(2.0)*kappa_pow2*spde.M1 + spde.M2;
}

/*
* Valid link functions
*/
enum valid_link {
  log_link                 = 0,
  logit_link               = 1,
  probit_link              = 2,
  inverse_link             = 3,
  identity_link            = 4,
  inverse_squared_link     = 5  // 1/x^2 default link with inverse gaussian
};
/*
* Valid families
*/
enum valid_family {
  gaussian_family = 0,
  binomial_family = 1,
  gamma_family    = 2,
  poisson_family  = 3,
  negative_binomial_family = 4
};
/*
* isNA
*/
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

/*
* Apply inverse link function for general cases
*/
template<class Type>
Type inverse_linkfun(Type& eta, int& link) {
  Type ans;
  switch (link) {
  case log_link:
    ans = exp(eta);
    break;
  case identity_link:
    ans = eta;
    break;
  case logit_link:
    ans = invlogit(eta);
    break;
  case probit_link:
    ans = pnorm(eta);
    break;
  case inverse_link:
    ans = Type(1) / eta;
    break;
  case inverse_squared_link:
    ans = Type(1) / (eta * eta);
    break;
    // TODO: Implement remaining links
  default:
    error("Link not implemented!");
  } // End switch
  return ans;
}

// transform Y -Inf-Inf -> X bound lb - ub
template <class Type> 
Type invlogit_general(Type Y, Type lb, Type ub) {
  return(lb + (ub - lb) * (1 / (1 + exp(-Y))));
}


template <class Type>
Type Gmean(vector<Type> x) {
  return(exp(log(x).sum() / x.size()));
}
