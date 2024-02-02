# @title update each effect once
# @param XtX a p by p matrix, X'X
# @param Xty a p vector
# @param s_init a list with elements sigma2, V, alpha, mu, Xr
# @param estimate_prior_variance boolean indicating whether to
#   estimate prior variance
# @param estimate_prior_method The method used for estimating prior
#   variance, 'optim' or 'EM'.
# @param check_null_threshold float a threshold on the log scale to
#   compare likelihood between current estimate and zero the null
# 
#' @importFrom Matrix diag
update_each_effect_ss = function (XtX, Xty, s_init,
                                  estimate_prior_variance = FALSE,
                                  estimate_prior_method = "optim",
                                  check_null_threshold = 0) {
  if (!estimate_prior_variance)
    estimate_prior_method = "none"
  
  # Repeat for each effect to update.
  s = s_init
  L = nrow(s$alpha)
  if (L > 0) {
    for (l in 1:L) {
      # Remove lth effect from fitted values.
      s$XtXr = s$XtXr - XtX %*% (s$alpha[l,] * s$mu[l,])
      # Compute residuals.
      XtR = Xty - s$XtXr
      # Correct zR discrepancy. For motivation see DENTIST paper Chen et al (2021)
      if (s$correct_zR_discrepancy$to_correct) {
        if (s$correct_zR_discrepancy$is_init) {
          # Skip check on the first iteration, but turn on the zR correction mode
          s$force_iterate = TRUE
          s$correct_zR_discrepancy$is_init = FALSE
        } else {
          # Get the current existing non-zero effect variables
          c_index = get_non_zero_effects_proxy(s$alpha, s$correct_zR_discrepancy$outlier_index) 
          if (length(c_index)>0) {
            # Detect outlier against existing non-zero effect variables
            outlier_index = detect_zR_discrepancy(c_index, s$correct_zR_discrepancy$outlier_index, 
                                                  XtR / sqrt(attr(XtX,"d")), XtX, r2=0.6, p=1E-4)
            # Apply correction
            if (outlier_index>0) {
              # cat(paste("\t- New outliers", paste(outlier_index), "for l =", l, "\n"))
              s$correct_zR_discrepancy$outlier_index = union(s$correct_zR_discrepancy$outlier_index, outlier_index)
              s$pi[s$correct_zR_discrepancy$outlier_index] = 0
              s$pi[s$pi>0]= s$pi[s$pi>0] / sum(s$pi)
            }
          }
        }
      }

      res = single_effect_regression_ss(as.matrix(XtR),attr(XtX,"d"),s$V[l],
              s$sigma2,s$pi,estimate_prior_method,check_null_threshold)
      
      # Update the variational estimate of the posterior mean.
      s$mu[l,]    = res$mu
      s$alpha[l,] = res$alpha
      s$mu2[l,]   = res$mu2
      s$V[l]      = res$V
      s$lbf[l]    = res$lbf_model
      s$lbf_variable[l,] = res$lbf
      s$KL[l]     = -res$lbf_model +
        SER_posterior_e_loglik_ss(attr(XtX,"d"),XtR,s$sigma2,
                                  res$alpha * res$mu,res$alpha * res$mu2)
      s$XtXr = s$XtXr + XtX %*% (s$alpha[l,] * s$mu[l,])
    }

    if (s$correct_zR_discrepancy$to_correct) {
      # Check if corrections are done, and if so let IBSS proceed as usual with convergence check
      # cat(paste("Removed mismatch at current iteration:", paste(sort(s$correct_zR_discrepancy$outlier_index), collapse=" "), "\n\n"))
      if (setequal(s$correct_zR_discrepancy$outlier_index,
                   s_init$correct_zR_discrepancy$outlier_index)) {
        s$correct_zR_discrepancy$outlier_stable_count = s$correct_zR_discrepancy$outlier_stable_count + 1
        if (s$correct_zR_discrepancy$outlier_stable_count >= s$correct_zR_discrepancy$outlier_stabilize) {
          # Wrap it up
          # cat("zR mismatch correction done!\n")
          s$correct_zR_discrepancy$to_correct = FALSE
          s$force_iterate = FALSE
        }
      } else {
        s$correct_zR_discrepancy$outlier_stable_count = 0 
      }
    }
  }
  s$XtXr = unname(as.matrix(s$XtXr))
  return(s)
}

detect_zR_discrepancy <- function(c_index, exclude_index, z, Rcov, r2=0.6, p=1E-4) {
  # > qchisq(1-1E-4,df=1)
  # [1] 15.13671
  # DENTIST-S test, $S(\hat{z}_1, \hat{z}_2, r_{12}) = \frac{(\hat{z}_1 - r_{12}\hat{z}_2)^2}{1-r_{12}^2} \sim \chi^2_{(1)}$
  dentist_s = function(z1,z2,r12) {
    (z1 - r12 * z2)^2 / (1 - r12^2)
  }
  is_sign_flip = function(z1,z2,r12) {
    ifelse(sign(z1) * sign(z2) * sign(r12) < 0, TRUE, FALSE)
  }

  # cat(paste("Non-zero set", paste(c_index, collapse=" "), "\n"))
  # every time here I want to just capture one outlier
  chisq_cutoff = qchisq(1-p, df = 1)
  fudge_factor = 1E-4
  x = abs(z)
  if (length(exclude_index)) {
    x[exclude_index] = -Inf
  }
  max_index = which.max(x)
  # cat(paste("\t- Test index:", max_index, "\n"))
  # FIXME: The 3 lines of codes commented out was an earlier version
  # where i assume that once a variant is determined causal it is not going to be considered outlier, in the current iteration
  # Although in later iterations it can still be removed, this might not be good enough because this bad variant can stay for quite a while.
  # However I've seen in the data that its immediately next buddy will revenge and push out the causal variant in the current iteration
  # Any since whatever I drop its not coming back, I may end up dropping anything that stands out after the fight, leaving its weaker buddies standing at the end.
  # which are not significant and lose power.
  # Still all considered I would rather be safe to keep FDR low as a priority. Therefore these 3 lines are commented out.
  #if (max_index %in% c_index) {
  #  return(-1)
  #}
  # Find the nearest correlation matrix from input
  # Because here our input is covariance
  # FIXME: is this correct?
  R = cov2cor(Rcov[c(max_index, c_index), c(max_index, c_index)]) * (1 - fudge_factor)
  diag(R) = 1
  z_test = z[c_index]
  z_max = z[max_index]
  stats_filter = sapply(1:length(z_test), function(i) dentist_s(z_max, z_test[i], R[1, i+1]))
  # cat(paste("\t- Chisq statistics", paste(round(stats_filter,3), collapse=" "), "\n"))
  stats_filter = (stats_filter > chisq_cutoff)
  r2_filter = sapply(1:length(z_test), function(i) R[1, i+1]^2)
  # cat(paste("\t- r2", paste(round(r2_filter,4), collapse=" "), "\n"))
  r2_filter = (r2_filter > r2)
  sign_filter = sapply(1:length(z_test), function(i) is_sign_flip(z_max, z_test[i], R[1, i+1]))
  # cat(paste("\t- Sign flip", paste(sign_filter, collapse=" "), "\n"))
  combined_filter = (stats_filter & (r2_filter | sign_filter))
  if(any(combined_filter)) {
    return (max_index)
  } else {
    return (-1)
  }
}

# KL Divergence
kl_divergence_p_q <- function(p, q) {
  p[which(p < .Machine$double.xmin)] <- .Machine$double.xmin
  q[which(q < .Machine$double.xmin)] <- .Machine$double.xmin
  sum(p * (log(p)-log(q)))
}

# Jensen-Shannon Divergence
js_divergence_p_q <- function(p, q) {
  m <- (p + q) / 2
  return((kl_divergence_p_q(p, m) + kl_divergence_p_q(q, m)) / 2)
}

get_non_zero_effects_proxy = function(alpha, exclude_index, tol=1E-4) {
  test_flat <- function(p0, tol) {
    p <- p0[p0 != 0]
    q <- 1 / length(p)
    js <- js_divergence_p_q(p, q)
    return(js < tol)
  }

  # effects to drop if alpha is flat
  to_drop <- which(apply(alpha, 1, test_flat, tol=tol))

  if (length(to_drop)) alpha = alpha[-to_drop,,drop=F]
  alpha [, exclude_index] = -Inf
  max_indices = apply(alpha, 1, FUN = which.max)
  return(unique(max_indices))
}