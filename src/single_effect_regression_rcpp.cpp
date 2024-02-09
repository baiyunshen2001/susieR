
#include <RcppArmadillo.h>
#include <iostream>
#include <stdexcept>
#ifdef _OPENMP
# include <omp.h>
#endif

using namespace Rcpp;
using namespace std;

using Rcpp::List;
using Rcpp::Named;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::IntegerVector;

using arma::vectorise;



NumericVector compute_tf_Xty_cpp(int order, const NumericVector& y) {
  NumericVector result = clone(y); // Make a copy of y to avoid modifying the original vector
  for (int i = 1; i <= order + 1; ++i) {
    NumericVector result_t=Rcpp::cumsum(result);
    for (int j = 1; j < result_t.size(); ++j)
      {
      result[j] = -1*result_t[j];
      }
  }
  return result;
}


NumericVector compute_Xty_cpp(const NumericMatrix& X, const NumericVector& y) {
  // Get attributes
  NumericVector cm = X.attr("scaled:center");
  NumericVector csd = X.attr("scaled:scale");

  // Convert to Armadillo structures
  arma::mat X_arma = as<arma::mat>(X);
  arma::vec y_arma = as<arma::vec>(y);

  // Compute cross product
  arma::vec ytX = arma::cross(y_arma, X_arma);

  // Check if X is a trend filtering matrix
  arma::vec scaled_Xty;
  if (!Rf_isNull(X.attr("matrix.type"))) {
    // Assuming compute_tf_Xty is defined
    int order = X.attr("order");
    scaled_Xty = compute_tf_Xty_cpp(order, y) / as<arma::vec>(csd);
  } else {
    scaled_Xty = ytX / as<arma::vec>(csd);
  }

  // Centering
  arma::vec centered_scaled_Xty = scaled_Xty - (as<arma::vec>(cm) / as<arma::vec>(csd)) * sum(y_arma);

  return wrap(centered_scaled_Xty);
}





double est_V_uniroot_cpp(NumericVector betahat, double shat2, NumericVector prior_weights) {
  // Define the R function 'uniroot'
  Function uniroot("uniroot");
  
  // Define the negloglik.grad.logscale function
  // Make sure this function is defined in R and available in the environment
  Function negloglik_grad_logscale("negloglik.grad.logscale");
  
  // Call the R 'uniroot' function
  // Note: 'negloglik_grad_logscale' is passed as an internal function
  List V_u = uniroot(Rcpp::_["f"] = negloglik_grad_logscale, 
                     Rcpp::_["interval"] = NumericVector::create(-10, 10),
                     Rcpp::_["extendInt"] = "upX",
                     Rcpp::_["betahat"] = betahat, 
                     Rcpp::_["shat2"] = shat2, 
                     Rcpp::_["prior_weights"] = prior_weights);
  
  // Extract the root and exponentiate
  double root = as<double>(V_u["root"]);
  return exp(root);
}



double neg_loglik_logscale(NumericVector V, NumericVector betahat, double shat2, NumericVector prior_weights) {
  
  int n = betahat.size();
  
  NumericVector lbf(n);
  
  // Calculate lbf for each SNP
  for(int i = 0; i < n; ++i) {
    if(std::isinf(shat2)) {
      lbf[i] = R::dnorm(betahat[i], 0, sqrt(V[i] + shat2), true) - R::dnorm(betahat[i], 0, sqrt(shat2), true);
    } else {
      lbf[i] = R::dnorm(betahat[i], 0, sqrt(V[i] + shat2), true) - R::dnorm(betahat[i], 0, sqrt(shat2), true);
    }
  }
  double maxlbf = max(lbf);
  NumericVector w = exp(lbf - maxlbf);
  NumericVector w_weighted = w * prior_weights;
  double weighted_sum_w = sum(w_weighted);
  
  // Return log likelihood
  return -(log(weighted_sum_w) + maxlbf);
}

double optimize_prior_variance_cpp(std::string optimize_V, NumericVector betahat, double shat2, 
                                   NumericVector prior_weights, Nullable<NumericVector> alpha, 
                                   Nullable<NumericVector> post_mean2, 
                                   Nullable<double> V_init, double check_null_threshold = 0){
  
  double V;
  if (V_init.isNotNull()) {
    V = as<double>(V_init);
  } else {
    V = 1.0; // Default value if V_init is NULL
  }
  
  std::string method = as<std::string>(optimize_V);
  
  if (method != "simple") {
    if (method == "optim") {
      // Define the R function
      Rcpp::Environment stats("package:stats"); 
      Rcpp::Function optim = stats["optim"];
      
      // Create initial value for optimization
      double initial_val = std::max(*std::max_element(betahat.begin(), betahat.end()) - shat2, 1.0);
      initial_val = log(initial_val);
      
      // Call the R optim function
      Rcpp::List opt_results = optim(Rcpp::_["par"] = initial_val,
                                     Rcpp::_["fn"] = Rcpp::InternalFunction(&neg_loglik_logscale),
                                     Rcpp::_["method"] = "Brent",
                                     Rcpp::_["lower"] = -30,
                                     Rcpp::_["upper"] = 15,
                                     Rcpp::_["betahat"] = betahat,
                                     Rcpp::_["shat2"] = shat2,
                                     Rcpp::_["prior_weights"] = prior_weights);
      
      // Extract the optimized parameter
      double lV = as<double>(opt_results["par"]);
      
      Function negloglik_grad_logscale("negloglik.grad.logscale");
      
      if(neg_loglik_logscale(Rcpp::_[lV]=lv,
                             Rcpp::_[betahat]=betahat, 
                             Rcpp::_[shat2]=shat2, 
                             Rcpp::_[prior_weights]=prior_weights) >
           neg_loglik_logscale(Rcpp::_[lV]=log(V), 
                               Rcpp::_[betahat]=betahat, 
                               Rcpp::_[shat2]=shat2, 
                               Rcpp::_[prior_weights]=prior_weights)){
        lV = log(V);
      }
      V = exp(lV);
    } else if (method == "uniroot") {
      V = est_V_uniroot_cpp(betahat, shat2, prior_weights);
    } else if (method == "EM") {
      if (alpha.isNotNull() && post_mean2.isNotNull()) {
        NumericVector alpha_vec = as<NumericVector>(alpha);
        NumericVector post_mean2_vec = as<NumericVector>(post_mean2);
        V = sum(alpha_vec * post_mean2_vec);
      }
    } else {
      stop("Invalid option for optimize_V method");
    }
  }
  
  return V;
}


// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List single_effect_regression_cpp(NumericVector y, NumericMatrix X, double V, 
                                  double residual_variance = 1, Nullable<NumericVector> prior_weights = R_NilValue,
                                  String optimize_V = "none", double check_null_threshold = 0) {
  
  
  NumericVector Xty = compute_Xty_cpp(X, y); // Placeholder for compute_Xty
  double d = X.attr("d");
  NumericVector betahat = (1 / d) * Xty;
  double shat2 = residual_variance / d;
  
  int ncolX = X.ncol();
  if (prior_weights.isNull()) {
    prior_weights = rep(1.0 / ncolX, ncolX);
  }
  
  
  if (optimize_V != "EM" && optimize_V != "none") {
    V = optimize_prior_variance_cpp(optimize_V, betahat, shat2, prior_weights, 
                                    NumericVector(), NumericVector(), V, 
                                    check_null_threshold);
  }
  Function dnorm("dnorm");
  
  NumericVector lbf(betahat.size());
  for (int i = 0; i < betahat.size(); ++i) {
    NumericVector lbf1 = dnorm(betahat[i], 0.0, sqrt(V + shat2), true);
    NumericVector lbf2 = dnorm(betahat[i], 0.0, sqrt(shat2), true);
    lbf[i] = lbf1[0] - lbf2[0];
  };
  
  lbf[std::isfinite(shat2)] = 0; // Handle special case of infinite shat2
  double maxlbf = max(lbf);
  NumericVector w = exp(lbf - maxlbf);
  NumericVector w_weighted = w * as<NumericVector>(prior_weights);
  double weighted_sum_w = sum(w_weighted);
  NumericVector alpha = w_weighted / weighted_sum_w;
  double post_var = 1.0 / (1.0 / V + d / residual_variance);
  NumericVector post_mean = (1.0 / residual_variance) * post_var * Xty;
  NumericVector post_mean2 = post_var + pow(post_mean, 2);
  
  double lbf_model = maxlbf + log(weighted_sum_w);
  double loglik = lbf_model + sum(dnorm(y, 0, sqrt(residual_variance), true));
  
  if (optimize_V == "EM") {
    V = optimize_prior_variance_cpp(optimize_V, betahat, shat2, prior_weights, 
                                    alpha, post_mean2, check_null_threshold);
  }
  
  return List::create(Named("alpha") = alpha, Named("mu") = post_mean,
                      Named("mu2") = post_mean2, Named("lbf") = lbf,
                      Named("lbf_model") = lbf_model, Named("V") = V,
                      Named("loglik") = loglik);
}

