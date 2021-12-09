/*
 * TMB source code for a univariate geostatistical CPUE model.
 * This source code applies nearest neighbour approach for the GF
 *
 */
#ifndef SpatialTemporalCPUE_NN_hpp
#define SpatialTemporalCPUE_NN_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

 /*
  *
  * The main obj$fn() call defines the joint negative log-likelihood
  */
template <class Type>
Type SpatialTemporalCPUENN(objective_function<Type>* obj) {
	// Dimensions
	DATA_INTEGER(n_i);         // Total number of observations
	DATA_INTEGER(n_t);         // Number of years
	/*
	 DATA_INTEGER( n_x );         // Number of vertices in SPDE mesh
	 DATA_INTEGER( n_z );         // Number of projection cells
	 DATA_INTEGER( p_s );         // number of spatial factors
	 DATA_INTEGER( p_c );         // number of catcspatialility factors
	 */
	 // Data
	DATA_VECTOR(y_i);       	 // data
	DATA_VECTOR(area_i);       	 // data
	DATA_IVECTOR(t_i);         // Time for each sample
	DATA_VECTOR(obs_t);        // observations per year
	DATA_IMATRIX(year_ndx_for_each_obs);   // C++ index, dim = [n_t, max(obs_t)]. -999 values indicate to stop reading, this will be a ragged matrix.
	DATA_VECTOR_INDICATOR(keep, y_i);       // NOTICE " keep " this is for OSA residuals
	DATA_STRUCT(spde, spde_t);            // spde object from INLA
	DATA_VECTOR(Proj_Area);              // Area for each projection grid
	DATA_INTEGER(family);                // 0 = gaussian, 1 = binomial, 2 = Gamma, 3 = Poisson see inst/HelperFuns.hpp for source code
	DATA_INTEGER(link);                  // 0 = log, 1 = logit, 2 = probit, 3 = inverse, 4 = identity
	DATA_MATRIX(time_model_matrix);      // year specific covariates [n_i, n_t]
	DATA_MATRIX(model_matrix);           // model matrix for non spatially related covariates (catcspatialility), month, vessel etc [n_i, p_c], intercept time_model_matrix
	// Model matrix timevariant spatial factors (Does not contain the intercept)
	DATA_ARRAY(X_spatial_ipt);	// Spatial projection model matrix for preference stuff dimensions = [n_i, p_s, n_t]
	DATA_ARRAY(X_spatial_proj_zpt);	// Spatial projection model matrix for preference stuff dimensions = [n_z, p_s, n_t]

	DATA_INTEGER(omega_indicator);       // is this component estimated in the model? 0 = no, 1 = yes
	DATA_INTEGER(epsilon_indicator);     // is this component estimated in the model? 0 = no, 1 = yes
	DATA_INTEGER(epsilon_ar1);             // 0 = no, 1 = yes, ignored if epsilon_indicator = 0
	DATA_IVECTOR(index_data_vertex);     // length n_i, links each data point to a vertex.
	DATA_IVECTOR(index_proj_vertex);     // length n_p, links each data point to a projecitons grid.

	// Preferential sampling
	DATA_VECTOR(pref_coef_bounds);  // pref_coef_bounds(0) = lower bound, pref_coef_bounds(1) = upper bound
	DATA_INTEGER(hyper_prior_on_pref);// 0 = no, 1 = yes
	DATA_INTEGER(apply_pref);       // 0 = no, 1 = yes
	DATA_ARRAY(Nij);                // number of observations in each each Proj cell dim[n_p, n_t]
	DATA_INTEGER(LCGP_approach);    // 0 Dinsdale approach, 1 = LGCP lattice approach
	// GAM stuff
	DATA_IVECTOR(spline_flag);                // length(2) spline_flag(0) = catcspatialility factors, spline_flag(1) spatial covariates, 0 no splines, 1 yes splines
	DATA_MATRIX(spline_model_matrix);        // model matrix catcspatialility coeffecients

	DATA_ARRAY(spline_spatial_model_matrix_ipt);        // model matrix catcspatialility coeffecients
	DATA_ARRAY(spline_spatial_model_matrix_proj_zpt);  // model matrix catcspatialility coeffecients

	DATA_SPARSE_MATRIX(S);                    //Penalization matrix diag(S1,S2,S3,S4,S5) without storing off-diagonal zeros.
	DATA_IVECTOR(Sdims);                      //Dimensions of S1,S2,S3,S4 and S5
	DATA_SPARSE_MATRIX(designMatrixForReport);//Design matrix for report of splines

	DATA_SPARSE_MATRIX(S_spatial);                    //Penalization matrix diag(S1,S2,S3,S4,S5) without storing off-diagonal zeros.
	DATA_IVECTOR(Sdims_spatial);                      //Dimensions of S1,S2,S3,S4 and S5
	DATA_SPARSE_MATRIX(designMatrixForReport_spatial);//Design matrix for report of splines

	DATA_IARRAY(spatial_constrained_coeff_ndx);       // rows are variables, columns with elements > 0 index coeffecients for spatial_constrained
	DATA_IVECTOR(spatial_covar_type);                 // one for each row of spatial_constrained_coeff_ndx, 0 = factor, 1 = numeric
	/* example structure
	 *  0  1   2   3        # categorical variable with 5 levels
	 *  4 -99 -99 -99       # numeric variable
	 *  5  6   7  -99       # categorical variable with 4 levels
	 */

	 // Simulate switches
	DATA_INTEGER(simulate_GF);                        // do you want to redraw a completely new GF? 0 = no, 1 = yes

	///////////////
	// Parameters
	///////////////
	// Beta fixed effects
	PARAMETER_VECTOR(betas);                          // p_c in length contains the intercept
	PARAMETER_VECTOR(constrained_spatial_betas);      // p_s - 1 in length the last coeffecient = -sum(spatial_betas)
	PARAMETER_VECTOR(constrained_time_betas);         // p_t = n_t in length coeffecients for time step

	// nuisance parameters
	PARAMETER(ln_phi);                                  // variance parameter for y_dist turn off for distributions such as Poisson.
	// Preference Model
	PARAMETER_VECTOR(logit_pref_coef);            // preferential parameter logistic transformation
	PARAMETER_VECTOR(logit_pref_hyper_params);  // vector continaing mu_pref and ln_sd_pref, in that order
	PARAMETER(lgcp_intercept);                          // intercept for the Log-Gaussian Cox Process Point process 

	// GMRF 
	PARAMETER(ln_kappa_omega);                  // spatial decay/range parameter if ln_kappa.size() = 1, then both omega and epsilon have the same range parameter
	PARAMETER(ln_tau_omega);                    // eta := kappa * tau, parameterisaton from https://github.com/nwfsc-assess/geostatistical_delta-GLMM/blob/master/inst/executables/geo_index_v4b.cpp
	PARAMETER(ln_kappa_epsilon);                // spatial decay/range parameter if ln_kappa.size() = 1, then both omega and epsilon have the same range parameter
	PARAMETER(ln_tau_epsilon);                  // eta := kappa * tau, parameterisaton from https://github.com/nwfsc-assess/geostatistical_delta-GLMM/blob/master/inst/executables/geo_index_v4b.cpp

	PARAMETER(trans_eps_rho);                     // Auto correlation parameter for Spatial temporal effect
	// Random effects - marginalised out
	PARAMETER_VECTOR(omega_input);              // spatial random effects for each vertix n_x
	PARAMETER_ARRAY(epsilon_input);             // spatial temporal random effects for each vertix [n_x x n_t]
	// GAM parameters
	PARAMETER_VECTOR(gammas);                           // Spline regression parameters
	PARAMETER_VECTOR(ln_lambda);                        // Penalization parameters
	PARAMETER_VECTOR(gammas_spatial);                       // Spline regression parameters for spatial coeffecients
	PARAMETER_VECTOR(ln_lambda_spatial);                    // Penalization parameters for spatial coeffecients

	///////////////////////////
	// Do some internal transformations
	///////////////////////////
	int i, t, k, j, obs_ndx;
	int n_p = Proj_Area.size();
	// un constrain spatial and catcspatialility coeffecients
	vector<Type> spatial_betas(X_spatial_ipt.col(0).cols());
	spatial_betas.setZero();
	Type constrained_totals = 0.0;
	int beta_ndx = 0;
	for (i = 0; i < spatial_constrained_coeff_ndx.rows(); ++i) {
		if (spatial_covar_type(i) == 1) {
			// numeric just a slope coeffecient no transformation needed
			spatial_betas(beta_ndx) = constrained_spatial_betas(spatial_constrained_coeff_ndx(i, 0));
			beta_ndx++;
		}
		else {
			constrained_totals = 0.0;
			for (j = 0; j < spatial_constrained_coeff_ndx.cols(); ++j) {
				// if 
				if ((j + 1) >= spatial_constrained_coeff_ndx.cols()) {
					spatial_betas(beta_ndx) = constrained_spatial_betas(spatial_constrained_coeff_ndx(i, j));
					constrained_totals += constrained_spatial_betas(spatial_constrained_coeff_ndx(i, j));
					++beta_ndx;
					spatial_betas(beta_ndx) = -1.0 * constrained_totals;
					++beta_ndx;
					break; // exit this covariate
				}
				else if (spatial_constrained_coeff_ndx(i, j) < 0) {
					spatial_betas(beta_ndx) = -1.0 * constrained_totals;
					++beta_ndx;
					break; // exit this covariate
				}
				else {
					spatial_betas(beta_ndx) = constrained_spatial_betas(spatial_constrained_coeff_ndx(i, j));
					constrained_totals += constrained_spatial_betas(spatial_constrained_coeff_ndx(i, j));
					++beta_ndx;
				}
			}
		}
	}

	vector<Type> time_betas(constrained_time_betas.size() + 1);
	for (i = 0; i < constrained_time_betas.size(); ++i)
		time_betas(i) = constrained_time_betas(i);
	time_betas(time_betas.size() - 1) = -1.0 * constrained_time_betas.sum();

	Type phi = exp(ln_phi);
	Type eps_rho = trans_eps_rho / sqrt(1.0 + pow(trans_eps_rho, 2));
	Type kappa_omega = exp(ln_kappa_omega);
	Type tau_omega = exp(ln_tau_omega);
	Type kappa_epsilon = exp(ln_kappa_epsilon);
	Type tau_epsilon = exp(ln_tau_epsilon);

	vector<Type> epsilon_vec(n_i);                                // temp vector
	vector<Type> spatial_Xbeta(n_i);
	vector<Type> spatial_splines(n_i);
	vector<Type> omega_i(n_i);//                                  // maps mesh to observations usign A
	vector<Type> epsilon_i(n_i);                                  // maps epsilon to observation
	vector<Type> spatial_covariate_i(n_i);                        // Spatail covariates linear
	vector<Type> spline_spatial_i(n_i);                                 // spatial spline contribution 
	vector<Type> pref_numerator(n_t);                             // sum^n_y (s(i))
	vector<Type> pref_denom(n_t);                                 // sum^N a_j exp(S(j))pref)
	vector<Type> relative_index(n_t);
	vector<Type> lambda = exp(ln_lambda);
	vector<Type> lambda_spatial = exp(ln_lambda_spatial);
	vector<Type> gamma_i;
	vector<Type> gamma_spatial_i;
	vector<Type> omega_proj(n_p);
	vector<Type> epsilon_proj(n_p);
	vector<Type> spatial_proj(n_p);
	SparseMatrix<Type> S_i;
	pref_numerator.setZero();
	pref_denom.setZero();
	spline_spatial_i.setZero();

	vector<Type> nll(7);  // 0 = GMRF (omega), 1 = GMRF (epsilon), 2 = obs, 3 = location, 4 = SPline catchability 5 = spline spatial, 6 = pref coeffecient hyper prior
	nll.setZero();


	// transform preferential sampling coeffecients
	vector<Type> pref_coef(n_t);
	if (logit_pref_coef.size() == 1) {
		for (t = 0; t < n_t; ++t)
			pref_coef(t) = invlogit_general(logit_pref_coef(0), pref_coef_bounds(0), pref_coef_bounds(1));
		if (hyper_prior_on_pref)
			nll(6) -= dnorm(logit_pref_coef(0), logit_pref_hyper_params(0), exp(logit_pref_hyper_params(1)));
	} else {
		for (t = 0; t < n_t; ++t) {
			pref_coef(t) = invlogit_general(logit_pref_coef(t), pref_coef_bounds(0), pref_coef_bounds(1));
			if (hyper_prior_on_pref)
				nll(6) -= dnorm(logit_pref_coef(t), logit_pref_hyper_params(0), exp(logit_pref_hyper_params(1)));
		}
	}

	// set some counters
	// Gaussian field stuff
	Type Range_omega = sqrt(8) / kappa_omega;
	Type MargSD_omega = 1 / sqrt(4 * M_PI) / tau_omega / kappa_omega;
	Type Range_epsilon = sqrt(8) / kappa_epsilon;
	Type MargSD_epsilon = 1 / sqrt(4 * M_PI) / tau_epsilon / kappa_epsilon;

	// Random field probability                                                                                                                              
	Eigen::SparseMatrix<Type> Q_omega = Q_spde_generalized(spde, kappa_omega, 2);
	Eigen::SparseMatrix<Type> Q_epsilon = Q_spde_generalized(spde, kappa_epsilon, 2);
	GMRF_t<Type> gmrf_omega = GMRF(Q_omega);
	GMRF_t<Type> gmrf_epsilon = GMRF(Q_epsilon);

	// GMRF same range for both omega and epsilon
	SparseMatrix<Type> Q = Q_spde(spde, kappa_omega);
	if (omega_indicator == 1) {
		SIMULATE{
		if (simulate_GF == 1) {
		omega_input = gmrf_omega.simulate();
		omega_input *= exp(-tau_omega);
		}
		}
		nll(0) = gmrf_omega(omega_input);                                 // Then the density is proportional to |Q|^.5*exp(-.5*x'*Q*x)
	}

	// Epsilon and spatial map to each observation 
	for (t = 0; t < epsilon_input.cols(); ++t) {
		if (epsilon_indicator == 1) {
			if (epsilon_ar1 == 1) {
				if (t == 0) {
					nll(1) += gmrf_epsilon(epsilon_input.col(t));
					SIMULATE{
						if (simulate_GF == 1) {
							epsilon_input.col(t) = gmrf_epsilon.simulate();
							epsilon_input.col(t) *= exp(-tau_epsilon);
						}
					}
				}
				else {
					nll(1) += gmrf_epsilon(epsilon_input.col(t) - eps_rho * epsilon_input.col(t - 1));
					SIMULATE{
						if (simulate_GF == 1) {
							epsilon_input.col(t) = gmrf_epsilon.simulate() + eps_rho * epsilon_input.col(t - 1);
							epsilon_input.col(t) *= exp(-tau_epsilon);
						}
					}
				}
			}
			else {
				nll(1) += gmrf_epsilon(epsilon_input.col(t));
				SIMULATE{
					if (simulate_GF == 1) {
					epsilon_input.col(t) = gmrf_epsilon.simulate();
					epsilon_input.col(t) *= exp(-tau_epsilon);
					}
				}
			}
		}
		epsilon_vec = vector<Type>(epsilon_input.col(t));
		spatial_Xbeta = X_spatial_ipt.col(t).matrix() * spatial_betas;
		spatial_splines = spline_spatial_model_matrix_ipt.col(t).matrix() * gammas_spatial;
		for (i = 0; i < obs_t(t); ++i) {
			obs_ndx = year_ndx_for_each_obs(t, i);
			spatial_covariate_i(obs_ndx) = spatial_Xbeta(obs_ndx);
			epsilon_i(obs_ndx) = epsilon_vec(index_data_vertex(obs_ndx)) / tau_epsilon;
			if (spline_flag(1) == 1)
				spline_spatial_i(obs_ndx) = spatial_splines(obs_ndx);
		}
	}


	// Spline stuff catcspatialility splines
	if (spline_flag(0) == 1) {
		k = 0;
		for (i = 0; i < Sdims.size(); i++) {
			int m_i = Sdims(i);
			gamma_i = gammas.segment(k, m_i);        // Recover betai
			S_i = S.block(k, k, m_i, m_i);      // Recover Si
			nll(4) -= Type(0.5) * m_i * ln_lambda(i) - 0.5 * lambda(i) * GMRF(S_i).Quadform(gamma_i);
			k += m_i;
		}
		vector<Type> splineForReport = designMatrixForReport * gammas;
		REPORT(splineForReport);
		ADREPORT(splineForReport);
		REPORT(lambda);
		REPORT(gammas);
	}
	// spatial splines
	if (spline_flag(1) == 1) {

		k = 0;
		for (i = 0; i < Sdims_spatial.size(); i++) {
			int m_i = Sdims_spatial(i);
			gamma_spatial_i = gammas_spatial.segment(k, m_i);        // Recover betai
			S_i = S_spatial.block(k, k, m_i, m_i);      // Recover Si
			nll(5) -= Type(0.5) * m_i * ln_lambda_spatial(i) - 0.5 * lambda_spatial(i) * GMRF(S_i).Quadform(gamma_spatial_i);
			k += m_i;
		}
		vector<Type> splineForReport_spatial = designMatrixForReport_spatial * gammas_spatial;
		REPORT(splineForReport_spatial);
		ADREPORT(splineForReport_spatial);
		REPORT(lambda_spatial);
		REPORT(gammas_spatial);
	}

	// Numerator for preference log likelihood
	for (i = 0; i < n_i; ++i) {
		omega_i(i) = omega_input(index_data_vertex(i)) / tau_omega;
		if (LCGP_approach == 0)
			pref_numerator(t_i(i)) += spatial_covariate_i(i) + spline_spatial_i(i) + epsilon_i(i) + omega_i(i);
	}
	// Systematic Component
	vector<Type> eta = model_matrix * betas + time_model_matrix * time_betas + spatial_covariate_i + spline_spatial_i + omega_i + epsilon_i;

	if (spline_flag(0) == 1)
		eta += spline_model_matrix * gammas;

	// Apply link function
	vector<Type> mu(eta.size());
	if (family == 1) {
		// Not going to area weight for binomial, doesn't make sense for proportion expected value
		for (i = 0; i < mu.size(); i++)
			mu(i) = inverse_linkfun(eta(i), link);
	}
	else {
		for (i = 0; i < mu.size(); i++)
			mu(i) = area_i(i) * inverse_linkfun(eta(i), link);
	}
	// Contribution Y | S, H
	// Observation likelihood
	Type s1, s2;
	Type tmp_loglik;
	for (i = 0; i < y_i.size(); i++) {
		if (!isNA(y_i(i))) {
			switch (family) {
			case gaussian_family:
				tmp_loglik = dnorm(y_i(i), mu(i), sqrt(phi), true);
				SIMULATE{ y_i(i) = rnorm(mu(i), sqrt(phi)); }
				break;
			case poisson_family:
				tmp_loglik = dpois(y_i(i), mu(i), true);
				SIMULATE{ y_i(i) = rpois(mu(i)); }
				break;
			case binomial_family:
				tmp_loglik = dbinom_robust(y_i(i), Type(1), eta(i), true); // Assumes logit link function Not probit
				SIMULATE{ y_i(i) = rbinom(Type(1), mu(i)); }
				break;
			case gamma_family:
				s1 = phi;           // shape
				s2 = mu(i) / phi;   // scale
				tmp_loglik = dgamma(y_i(i), s1, s2, true);
				SIMULATE{ y_i(i) = rgamma(s1, s2); }
				break;
			case negative_binomial_family:
				s1 = log(mu(i));                       // log(mu)
				s2 = 2. * s1 - ln_phi;                         // log(var - mu)
				tmp_loglik = dnbinom_robust(y_i(i), s1, s2, true);
				SIMULATE{
				s1 = mu(i);
				s2 = mu(i) * (1.0 + phi);  // (1+phi) guarantees that var >= mu
				y_i(i) = rnbinom2(s1, s2);
				}
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
	for (i = 0; i < n_p; ++i)
		omega_proj(i) = omega_input(index_proj_vertex(i)) / tau_omega;

	// time-varying components

	for (t = 0; t < n_t; ++t) {
		for (i = 0; i < n_p; ++i)
			epsilon_proj(i) = epsilon_input(index_proj_vertex(i), t) / tau_epsilon;

		spatial_proj = X_spatial_proj_zpt.col(t).matrix() * spatial_betas + omega_proj + epsilon_proj;
		if (spline_flag(1) == 1)
			spatial_proj += spline_spatial_model_matrix_proj_zpt.col(t).matrix() * gammas_spatial;

		relative_index(t) = (Proj_Area * exp(time_betas(t) + epsilon_proj)).sum();

		if (LCGP_approach == 0)
			pref_denom(t) = log((Proj_Area * exp(pref_coef(t) * spatial_proj)).sum());
		// deal with Pr(S | W, X) Observations
		if (apply_pref == 1) {
			if (LCGP_approach == 0) {
				nll(3) -= pref_coef(t) * pref_numerator(t) - obs_t(t) * pref_denom(t);
			}
			else if (LCGP_approach == 1) {
				for (i = 0; i < n_p; ++i)
					nll(3) -= dpois(Nij(i, t), Proj_Area(i) * exp(lgcp_intercept + pref_coef(t) * spatial_proj(i)), true);
			}
		}
	}


	// Time-serie
	Type gmean = Gmean(relative_index);
	vector<Type> standardised_index = relative_index / Gmean(relative_index);
	//for(i = 0; i < nll.size(); ++i)
	//std::cout << nll(i) << " ";
	//std::cout << std::endl;

	SIMULATE{
	REPORT(y_i);
	}
	vector<Type> betas_w_intercept(betas.size());
	betas_w_intercept(0) = betas(0);
	for (i = 1; i < betas_w_intercept.size(); ++i)
		betas_w_intercept(i) = betas(0) + betas(i);

	//////////////////
	// Report section
	//////////////////
	REPORT(MargSD_omega);
	REPORT(MargSD_epsilon);

	REPORT(Range_omega);
	REPORT(Range_epsilon);

	REPORT(tau_omega);
	REPORT(tau_epsilon);

	REPORT(kappa_omega);
	REPORT(kappa_epsilon);

	REPORT(pref_numerator);
	REPORT(pref_denom);

	REPORT(nll);
	REPORT(betas);
	REPORT(spatial_betas);
	REPORT(spatial_splines)
		REPORT(time_betas);
	REPORT(betas_w_intercept);
	REPORT(relative_index);
	//REPORT( standardised_index );
	//REPORT( gmean );

	REPORT(omega_input);
	REPORT(epsilon_input);
	REPORT(omega_i);
	REPORT(epsilon_i);
	REPORT(spatial_covariate_i);
	REPORT(spline_spatial_i);

	REPORT(omega_proj);
	//REPORT( epsilon_proj ); // shouldn't report as these are time dependent so wouldn't make sense to report this to users
	//REPORT( spatial_proj );

	REPORT(pref_coef);
	REPORT(lgcp_intercept);
	REPORT(logit_pref_hyper_params);

	REPORT(phi);
	REPORT(mu);
	REPORT(eta);
	REPORT(eps_rho);

	// ADREPORT
	ADREPORT(phi);
	ADREPORT(pref_coef);
	ADREPORT(Range_omega);
	ADREPORT(Range_epsilon);
	ADREPORT(MargSD_omega);
	ADREPORT(MargSD_epsilon);
	ADREPORT(spatial_betas);
	ADREPORT(betas);
	ADREPORT(time_betas);
	ADREPORT(betas_w_intercept);

	ADREPORT(relative_index);
	ADREPORT(standardised_index);
	ADREPORT(eps_rho);

	return nll.sum();
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif