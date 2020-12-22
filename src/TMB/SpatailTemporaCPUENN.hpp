/// @file ModelA.hpp

// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

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

// name of function below **MUST** match filename
// (here it's ModelA)
template <class Type>
Type SpatialTemporalCPUENN(objective_function<Type>* obj) {
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  // Dimensions
  DATA_INTEGER( n_i );         // Total number of observations
  DATA_INTEGER( n_t );         // Number of years
  DATA_INTEGER( n_p );         // Number of projection grids
  
  /*
   DATA_INTEGER( n_x );         // Number of vertices in SPDE mesh
   DATA_INTEGER( n_z );         // Number of projection cells
   DATA_INTEGER( p_s );         // number of spatial factors
   DATA_INTEGER( p_c );         // number of catchability factors
   */
  // Data
  DATA_VECTOR( y_i );       	 // data
  DATA_VECTOR( area_i );       	 // data
  DATA_IVECTOR( t_i );         // Time for each sample
  DATA_VECTOR( obs_t );        // observations per year
  DATA_IMATRIX( year_ndx_for_each_obs );   // C++ index, dim = [n_t, max(obs_t)]. -999 values indicate to stop reading, this will be a ragged matrix.
  DATA_VECTOR_INDICATOR (keep , y_i);       // NOTICE " keep " this is for OSA residuals
  //DATA_SPARSE_MATRIX( A ) ;               // Sparse matrix that maps, nodes to each observation location
  DATA_STRUCT( spde, spde_t );            // spde object from INLA
  DATA_VECTOR( Proj_Area ) ;              // Area for each projection grid
  DATA_INTEGER( family ) ;                // 0 = gaussian, 1 = binomial, 2 = Gamma, 3 = Poisson
  DATA_INTEGER( link ) ;                  // 0 = log, 1 = logit, 2 = probit, 3 = inverse, 4 = identity
  DATA_MATRIX( time_model_matrix ) ;      // year specific covariates [n_i, n_t]
  DATA_MATRIX( model_matrix ) ;           // model matrix for non spatially related covariates (catchability), month, vessel etc [n_i, p_c], intercept time_model_matrix
  // Model matrix timevariant spatial factors (Does not contain the intercept)
  DATA_ARRAY( X_spatial_ipt );	// Spatial projection model matrix for preference stuff dimensions = [n_i, p_s, n_t]
  DATA_ARRAY( X_spatial_proj_zpt );	// Spatial projection model matrix for preference stuff dimensions = [n_z, p_s, n_t]
  
  DATA_VECTOR( pref_coef_bounds );  // pref_coef_bounds(0) = lower bound, pref_coef_bounds(1) = upper bound
  DATA_INTEGER( apply_pref );       // 0 = no, 1 = yes
  
  DATA_INTEGER( omega_indicator );       // 0 = no, 1 = yes
  DATA_INTEGER( epsilon_indicator );     // 0 = no, 1 = yes
  DATA_IVECTOR( index_data_vertex );     // length n_i, links each data point to a vertex.
  DATA_IVECTOR( index_proj_vertex );     // length n_p, links each data point to a projecitons grid.
  
  // GAM stuff
  DATA_IVECTOR( spline_flag );                // length(2) spline_flag(0) = catchability factors, spline_flag(1) habitat covariates, 0 no splines, 1 yes splines
  DATA_MATRIX( spline_model_matrix ) ;        // model matrix catchability coeffecients
  
  DATA_ARRAY( spline_habitat_model_matrix_ipt ) ;        // model matrix catchability coeffecients
  DATA_ARRAY( spline_habitat_model_matrix_proj_zpt ) ;  // model matrix catchability coeffecients
  
  DATA_SPARSE_MATRIX( S );                    //Penalization matrix diag(S1,S2,S3,S4,S5) without storing off-diagonal zeros.
  DATA_IVECTOR( Sdims );                      //Dimensions of S1,S2,S3,S4 and S5
  DATA_SPARSE_MATRIX( designMatrixForReport );//Design matrix for report of splines
  
  DATA_SPARSE_MATRIX( S_hab );                    //Penalization matrix diag(S1,S2,S3,S4,S5) without storing off-diagonal zeros.
  DATA_IVECTOR( Sdims_hab );                      //Dimensions of S1,S2,S3,S4 and S5
  DATA_SPARSE_MATRIX( designMatrixForReport_hab );//Design matrix for report of splines
  
  
  ///////////////
  // Parameters
  ///////////////
  // Beta fixed effects
  PARAMETER_VECTOR( betas );                          // p_c in length contains the intercept
  PARAMETER_VECTOR( constrained_spatial_betas );      // p_s - 1 in length the last coeffecient = -sum(spatial_betas)
  PARAMETER_VECTOR( constrained_time_betas );         // p_t = n_t in length coeffecients for time step
  
  // nuisance parameters
  PARAMETER(ln_phi);                      // variance parameter for y_dist turn off for distributions such as Poisson.
  // Preference Model
  PARAMETER(logit_pref_coef);                    // preferential parameter
  // GMRF 
  PARAMETER( ln_kappa );                  // spatial decay/range parameter if ln_kappa.size() = 1, then both omega and epsilon have the same range parameter
  PARAMETER( ln_eta_omega );              // eta := kappa * tau, parameterisaton from https://github.com/nwfsc-assess/geostatistical_delta-GLMM/blob/master/inst/executables/geo_index_v4b.cpp
  PARAMETER( ln_eta_epsilon );            // eta := kappa * tau
  
  PARAMETER_VECTOR( ln_kappa );           // spatial decay/range parameter if ln_kappa.size() = 1, then both omega and epsilon have the same range parameter
  PARAMETER_VECTOR( ln_tau );           // tau
  
  
  
  PARAMETER(logit_eps_rho);               // Auto correlation parameter for Spatial temporal effect
  // Random effects - marginalised out
  PARAMETER_VECTOR( omega_input );              // spatial random effects for each vertix n_x
  PARAMETER_ARRAY( epsilon_input );             // spatial temporal random effects for each vertix [n_x x n_t]
  // GAM parameters
  PARAMETER_VECTOR( gammas );                           // Spline regression parameters
  PARAMETER_VECTOR( ln_lambda );                        // Penalization parameters
  PARAMETER_VECTOR( gammas_hab );                       // Spline regression parameters for habitat coeffecients
  PARAMETER_VECTOR( ln_lambda_hab );                    // Penalization parameters for habitat coeffecients
  
  ///////////////////////////
  // Do some internal transformations
  ///////////////////////////
  // un constrain spatial and catchability coeffecients
  vector<Type> spatial_betas(constrained_spatial_betas.size() + 1);
  vector<Type> time_betas(constrained_time_betas.size() + 1);
  
  for(int i = 0; i < constrained_spatial_betas.size(); ++i) 
    spatial_betas(i) = constrained_spatial_betas(i);
  spatial_betas(spatial_betas.size() - 1) = -1.0 * constrained_spatial_betas.sum();
  for(int i = 0; i < constrained_time_betas.size(); ++i) 
    time_betas(i) = constrained_time_betas(i);
  time_betas(time_betas.size() - 1) = -1.0 * constrained_time_betas.sum();
  
  Type phi = exp(ln_phi); 
  Type ln_tau_omega = ln_eta_omega - ln_kappa;
  Type ln_tau_epsilon = ln_eta_epsilon - ln_kappa;
  Type tau_omega = exp(ln_tau_omega);
  Type tau_epsilon = exp(ln_tau_epsilon);
  Type kappa = exp(ln_kappa);             // Range
  Type pref_coef = invlogit_general(logit_pref_coef, pref_coef_bounds(0), pref_coef_bounds(1));
  Type eps_rho = invlogit_general(logit_eps_rho, Type(-0.99), Type(0.99));
  vector<Type> epsilon_vec(n_i);                                // temp vector
  vector<Type> spatial_Xbeta(n_i);
  vector<Type> spatial_splines(n_i);
  vector<Type> omega_i(n_i);//                     // maps mesh to observations usign A
  vector<Type> epsilon_i(n_i);                                  // maps epsilon to observation
  vector<Type> spatial_covariate_i(n_i);                        // Spatail covariates linear
  vector<Type> spline_hab_i(n_i);                                 // Habitat spline contribution 
  vector<Type> pref_numerator(n_t);                             // sum^n_y (s(i))
  vector<Type> pref_denom(n_t);                                 // sum^N a_j exp(S(j))pref)
  vector<Type> relative_index(n_t);
  vector<Type> lambda = exp(ln_lambda);
  vector<Type> lambda_hab = exp(ln_lambda_hab);
  vector<Type> gamma_i(n_i);
  vector<Type> gamma_hab_i(n_i);
  vector<Type> omega_proj(n_p);
  vector<Type> eta(n_i);
  vector<Type> mu(n_i);
  vector<Type> epsilon_proj(n_p);
  vector<Type> spatial_proj(n_p);
  SparseMatrix<Type> S_i; 
  pref_numerator.setZero();
  pref_denom.setZero();
  spline_hab_i.setZero();
  
  vector<Type> nll(6);  // 0 = GMRF (omega), 1 = GMRF (epsilon), 2 = obs, 3 = location, 4 = SPline catchability 5 = spline habitat
  nll.setZero();
  // set some counters
  int i, t, k, obs_ndx;
  
  // GMRF same range for both omega and epsilon
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  GMRF_t<Type> gmrf_Q = GMRF(Q);
  // Omega nll contribution
  if (omega_indicator == 1)
    nll(0) = SCALE(GMRF(Q),  Type(1.0)/ tau_omega)(omega_input);                                 // Then the density is proportional to |Q|^.5*exp(-.5*x'*Q*x)
  
  // Epsilon and spatial map to each observation 
  for(t = 0; t < epsilon_input.cols(); ++t) {
    epsilon_vec = vector<Type>(epsilon_input.col(t));
    spatial_Xbeta = X_spatial_ipt.col(t).matrix() * spatial_betas;
    spatial_splines = spline_habitat_model_matrix_ipt.col(t).matrix() * gammas_hab;
    for(i = 0; i < obs_t(t); ++i) {
      obs_ndx = year_ndx_for_each_obs(t,i);
      spatial_covariate_i(obs_ndx) = spatial_Xbeta(obs_ndx);
      epsilon_i(obs_ndx) = epsilon_vec(index_data_vertex(obs_ndx));
      if (spline_flag(1) == 1) 
        spline_hab_i(obs_ndx) = spatial_splines(obs_ndx);
    }
    
    // Omega nll contribution
    if (epsilon_indicator == 1)
      nll(1) += SCALE(gmrf_Q, Type(1.0) / tau_epsilon)(epsilon_input.col(t));
    //} else {
    //  nll(1) += gmrf_Q(epsilon_input.col(t) - eps_rho * epsilon_input.col(t-1));
    //}
  }
   
  //nll(1) = SEPARABLE(AR1(eps_rho), GMRF(Q))(Epsilon_input);   // AR(1)
  
  // Spline stuff catchability splines
  
  if (spline_flag(0) == 1) {
    k = 0;
    for( i = 0;i < Sdims.size(); i++){
      int m_i = Sdims(i);
      gamma_i = gammas.segment(k,m_i);        // Recover betai
      S_i = S.block(k,k,m_i,m_i);      // Recover Si
      nll(4) -= Type(0.5) * m_i * ln_lambda(i) - 0.5 * lambda(i) * GMRF(S_i).Quadform(gamma_i);
      k += m_i;
    }
    vector<Type> splineForReport = designMatrixForReport * gammas;
    REPORT( splineForReport );
    REPORT( lambda );
    REPORT( gammas );
  }
  // Habitat splines
  if (spline_flag(1) == 1) {
    
    k = 0;
    for( i = 0;i < Sdims_hab.size(); i++){
      int m_i = Sdims_hab(i);
      gamma_hab_i = gammas_hab.segment(k,m_i);        // Recover betai
      S_i = S_hab.block(k,k,m_i,m_i);      // Recover Si
      nll(5) -= Type(0.5) * m_i * ln_lambda_hab(i) - 0.5 * lambda_hab(i) * GMRF(S_i).Quadform(gamma_hab_i);
      k += m_i;
    }
    vector<Type> splineForReport_hab = designMatrixForReport_hab * gammas_hab;
    REPORT( splineForReport_hab );
    REPORT( lambda_hab );
    REPORT( gammas_hab );
  }
  
  // Numerator for preference log likelihood
  for(i = 0; i < n_i; ++i) {
    omega_i(i) = omega_input(index_data_vertex(i));
    pref_numerator(t_i(i)) += spatial_covariate_i(i) + spline_hab_i(i) + epsilon_i(i) + omega_i(i);
  }
  // Systematic Component
  eta =  model_matrix * betas + time_model_matrix * time_betas + spatial_covariate_i + spline_hab_i + omega_i + epsilon_i;
  
  if (spline_flag(0) == 1)
    eta += spline_model_matrix * gammas;
  
  // Apply link function
  //vector<Type> mu(eta.size());
  for (i = 0; i < mu.size(); i++)
    mu(i) = area_i(i) * inverse_linkfun(eta(i), link);
  
  // Contribution Y | S, H
  // Observation likelihood
  Type s1, s2;
  Type tmp_loglik;
  for (i = 0; i < y_i.size(); i++){
    if ( !isNA(y_i(i)) ) {
      switch (family) {
      case gaussian_family:
        tmp_loglik = dnorm(y_i(i), mu(i), sqrt(phi), true);
        SIMULATE{y_i(i) = rnorm(mu(i), sqrt(phi));}
        break;
      case poisson_family:
        tmp_loglik = dpois(y_i(i), mu(i), true);
        SIMULATE{y_i(i) = rpois(mu(i));}
        break;
      case binomial_family:
        tmp_loglik = dbinom_robust(y_i(i), Type(1), mu(i), true);
        SIMULATE{y_i(i) = rbinom(Type(1), mu(i));}
        break;
      case Gamma_family:
        s1 = phi;           // shape
        s2 = mu(i) / phi;   // scale
        tmp_loglik = dgamma(y_i(i), s1, s2, true);
        SIMULATE{y_i(i) = rgamma(s1, s2);}
        break;
      default:
        error("Family not implemented!");
      } // End switch
      // Add up
      nll(2) -= keep(i) * tmp_loglik;
    }
  }
  
  
  // deal with location density and projections
  // Location density
  // Projections
  for(i = 0; i < n_p; ++i) {
    omega_proj(i) = omega_input(index_proj_vertex(i));
  }
  
  // time-varying components
  for(t = 0; t < n_t; ++t) {
    for(i = 0;i < n_p; ++i)
      epsilon_proj(i) = epsilon_input(index_proj_vertex(i), t);
    
    spatial_proj = X_spatial_proj_zpt.col(t).matrix() * spatial_betas + omega_proj + epsilon_proj;
    if (spline_flag(1) == 1)
      spatial_proj += spline_habitat_model_matrix_proj_zpt.col(t).matrix() * gammas_hab;
    
    relative_index(t) = (Proj_Area * exp(time_betas(t) + spatial_proj)).sum();
    
    
    pref_denom(t) = log((Proj_Area * exp(pref_coef * spatial_proj)).sum()); 
    // deal with Pr(S | W, X) Observations
    if(apply_pref == 1) {
      nll(3) -=  pref_coef * pref_numerator(t) - obs_t(t) * pref_denom(t);
    }
  }
  
  
  // GMRF variables decorrelation and marginal std deviation
  Type Range = sqrt(8) / kappa;
  Type SD_omega = 1 / sqrt(4 * M_PI * exp(2*ln_tau_omega) * exp(2 * ln_kappa));
  Type SD_epsilon = 1 / sqrt(4 * M_PI * exp(2*ln_tau_epsilon) * exp(2 * ln_kappa));
  
  // Time-serie
  Type gmean = Gmean(relative_index);
  vector<Type> standardised_index = relative_index / Gmean(relative_index);
  
  /*
  for(i = 0; i < nll.size(); ++i)
  std::cout << nll(i) << " ";
  std::cout << std::endl;
  */
  
  //////////////////
  // Report section
  //////////////////
  REPORT( SD_epsilon );
  REPORT( SD_omega );
  REPORT( Range );
  REPORT( tau_omega );
  REPORT( tau_epsilon );
  REPORT( kappa );
  
  
  REPORT( pref_numerator );
  REPORT( pref_denom );
  
  REPORT( nll );
  REPORT( betas );
  REPORT( spatial_betas );
  REPORT( spatial_splines )
  REPORT( time_betas );
  REPORT( relative_index);
  //REPORT( standardised_index );
  //REPORT( gmean );
  
  REPORT( omega_input );
  REPORT( epsilon_input );
  REPORT( omega_i );
  REPORT( epsilon_i );
  REPORT( spatial_covariate_i );
  REPORT( spline_hab_i );
  REPORT( omega_proj );
  REPORT( epsilon_proj );
  REPORT( spatial_proj );
  
  REPORT( pref_coef );
  
  REPORT( phi );
  REPORT( mu );
  REPORT( eta );
  // ADREPORT
  ADREPORT( phi );
  ADREPORT( Range );
  ADREPORT( pref_coef );
  ADREPORT( SD_epsilon );
  ADREPORT( SD_omega );
  //ADREPORT( relative_index );
  ADREPORT( standardised_index );
  
  return nll.sum();
}

#undef TMB_OBJECTIVE_PTR

