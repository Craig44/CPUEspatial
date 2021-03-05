/*
* Valid link functions
*/
enum valid_link {
  log_link                 = 0,
  logit_link               = 1,
  probit_link              = 2,
  inverse_link             = 3,
  identity_link            = 4
};
/*
* Valid families
*/
enum valid_family {
  gaussian_family = 0,
  binomial_family = 1,
  Gamma_family    = 2,
  poisson_family  = 3
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
