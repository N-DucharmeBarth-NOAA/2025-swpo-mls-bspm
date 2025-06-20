quick_diagnostics <- function(fit) {
  cat("=== STAN FIT DIAGNOSTICS ===\n")
  
  # Basic info
  cat("Chains:", fit@sim$chains, "\n")
  cat("Iterations:", fit@sim$iter, "\n")
  cat("Warmup:", fit@sim$warmup, "\n")
  
  # Convergence
  rhats <- summary(fit)$summary[,"Rhat"]
  cat("Max Rhat:", round(max(rhats, na.rm = TRUE), 4), "\n")
  cat("Parameters with Rhat > 1.01:", sum(rhats > 1.01, na.rm = TRUE), "\n")
  
  # Effective sample size
  n_eff <- summary(fit)$summary[,"n_eff"]
  cat("Min n_eff:", round(min(n_eff, na.rm = TRUE), 0), "\n")
  
  # Divergences
  div <- sum(get_divergent_iterations(fit))
  cat("Divergent transitions:", div, "\n")
  
  # Energy diagnostics
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  energy <- sapply(sampler_params, function(x) x[,"energy__"])
  cat("Energy diagnostics check:\n")
  check_energy(fit)
  
  cat("========================\n")
}

compare_marginals <- function(fit1, fit2, params = c("logK", "r", "shape", "sigmap", "sigmao_add", "qeff", "rho", "sigma_qdev")) {
  for(param in params) {
    samples1 <- extract(fit1, param)[[1]]
    samples2 <- extract(fit2, param)[[1]]
    
    cat("\n", param, ":\n")
    cat("  Model 1 (ppc) - Mean:", round(mean(samples1), 3), 
        "SD:", round(sd(samples1), 3), "\n")
    cat("  Model 2 (fit) - Mean:", round(mean(samples2), 3), 
        "SD:", round(sd(samples2), 3), "\n")
    cat("  Difference in means:", round(mean(samples2) - mean(samples1), 3), "\n")
  }
}
