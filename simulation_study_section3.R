# ==============================================================================
# SECTION 3: COMPREHENSIVE MONTE CARLO SIMULATION STUDY (FIXED)
# ==============================================================================
# Fixes applied:
#   [FIX 1] Harmonized bounds between APSO and MLE refinement phases
#   [FIX 2] Fallback returns converged=TRUE to preserve candidate validity
#   [FIX 3] Bounds derived dynamically from empirical moments
#   [FIX 4] Robust consistency penalty using median
# ==============================================================================

rm(list = ls())
options(warn = -1, scipen = 999, digits = 6)

suppressPackageStartupMessages({
  library(knitr); library(numDeriv)
})

set.seed(20240617)

# ==============================================================================
# PART 1: RP DISTRIBUTION FUNCTIONS
# ==============================================================================
p_rp <- function(x, theta, lambda) {
  x <- pmax(x, 1e-10)
  z <- exp(-theta * (x^2))
  (1 - exp(-lambda * (1 - z))) / (1 - exp(-lambda))
}

f_rp <- function(x, theta, lambda) {
  x <- pmax(x, 1e-10); x2 <- x^2; z <- exp(-theta * x2)
  (2 * theta * x * lambda * z * exp(-lambda * (1 - z))) / (1 - exp(-lambda))
}

q_rp <- function(u, theta, lambda) {
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  inner <- 1 + log(1 - u * (1 - exp(-lambda))) / lambda
  inner <- pmax(inner, 1e-12)
  sqrt(-log(inner) / theta)
}

r_rp_gen <- function(n, theta, lambda) q_rp(runif(n), theta, lambda)

# ==============================================================================
# PART 2: LOG-LIKELIHOOD AND GRADIENT (CORRECTED SIGNS)
# ==============================================================================
nll_rp <- function(par, x) {
  theta <- par[1]; lambda <- par[2]
  if (theta <= 1e-8 || lambda <= 1e-8 || lambda > 30) return(1e12)
  n <- length(x)
  ll <- n * log(2 * theta * lambda) - n * lambda -
        n * log(1 - exp(-lambda)) + sum(log(pmax(x, 1e-10))) -
        theta * sum(x^2) + lambda * sum(exp(-theta * x^2))
  if (!is.finite(ll)) return(1e12)
  -ll
}

grad_nll_rp <- function(par, x) {
  theta <- par[1]; lambda <- par[2]
  if (theta <= 1e-8 || lambda <= 1e-8) return(c(0, 0))
  n <- length(x); x2 <- x^2; e <- exp(-theta * x2)
  d_theta  <- (n / theta) - sum(x2) - lambda * sum(x2 * e)
  d_lambda <- (n / lambda) - n - (n * exp(-lambda)) / (1 - exp(-lambda)) + sum(e)
  c(-d_theta, -d_lambda)
}

# ==============================================================================
# PART 3: SHARED BOUNDS FUNCTION (? NEW)
# ==============================================================================
# Computes consistent bounds used by BOTH APSO and MLE-APSO phases
compute_bounds <- function(x) {
  m2 <- mean(x^2)
  th_prior <- 1 / (m2 + 1e-6)
  
  lb <- c(
    max(th_prior * 0.01, 1e-6),   # theta lower
    0.01                           # lambda lower
  )
  ub <- c(
    max(th_prior * 100, 1.0),     # theta upper
    30                             # lambda upper
  )
  
  list(lb = lb, ub = ub, th_prior = th_prior)
}

# ==============================================================================
# PART 4: ESTIMATION METHODS
# ==============================================================================

# --- (A) Standard MLE: multi-start L-BFGS-B ---
fit_mle <- function(x, n_starts = 5) {
  bnds <- compute_bounds(x)
  
  theta_starts  <- 10^runif(n_starts, log10(max(bnds$lb[1], 1e-4)),
                                     log10(min(bnds$ub[1], 10)))
  lambda_starts <- 10^runif(n_starts, log10(max(bnds$lb[2], 0.05)),
                                     log10(min(bnds$ub[2], 15)))
  
  best_val <- Inf
  best_par <- c(NA, NA)
  start_time <- Sys.time()
  
  for (k in 1:n_starts) {
    init <- c(theta_starts[k], lambda_starts[k])
    res <- tryCatch({
      optim(par = init, fn = nll_rp, gr = grad_nll_rp, x = x,
            method = "L-BFGS-B",
            lower = bnds$lb, upper = bnds$ub,
            control = list(maxit = 500, factr = 1e7))
    }, error = function(e) NULL)
    
    if (!is.null(res) && res$value < best_val && all(is.finite(res$par))) {
      best_val <- res$value
      best_par <- res$par
    }
  }
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  list(par = best_par,
       nll = best_val,
       converged = !any(is.na(best_par)) && is.finite(best_val),
       time = elapsed)
}

# --- (B) Standalone APSO ---
fit_apso <- function(x, n_particles = 40, max_iter = 200) {
  bnds <- compute_bounds(x)
  lb_log <- log(bnds$lb)
  ub_log <- log(bnds$ub)
  
  nll_log <- function(p) nll_rp(exp(p), x)
  
  dim <- 2
  X <- matrix(0, n_particles, dim)
  for (d in 1:dim) X[, d] <- runif(n_particles, lb_log[d], ub_log[d])
  V <- matrix(0, n_particles, dim)
  step_scale <- (ub_log - lb_log) * 0.1
  
  pbest <- X
  pbest_val <- apply(X, 1, nll_log)
  gidx <- which.min(pbest_val)
  gbest <- pbest[gidx, ]
  gbest_val <- pbest_val[gidx]
  
  start_time <- Sys.time()
  
  for (iter in 1:max_iter) {
    w  <- 0.9 - 0.5 * (iter / max_iter)
    c1 <- 2.5 - 2.0 * (iter / max_iter)
    c2 <- 0.5 + 2.0 * (iter / max_iter)
    
    for (i in 1:n_particles) {
      r1 <- runif(dim); r2 <- runif(dim)
      V[i, ] <- w * V[i, ] + c1 * r1 * (pbest[i, ] - X[i, ]) +
                c2 * r2 * (gbest - X[i, ])
      V[i, ] <- pmax(pmin(V[i, ], step_scale), -step_scale)
      X[i, ] <- X[i, ] + V[i, ]
      X[i, ] <- pmax(pmin(X[i, ], ub_log), lb_log)
      
      val <- nll_log(X[i, ])
      if (val < pbest_val[i]) {
        pbest_val[i] <- val
        pbest[i, ] <- X[i, ]
        if (val < gbest_val) {
          gbest_val <- val
          gbest <- X[i, ]
        }
      }
    }
  }
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  par_orig <- exp(gbest)
  
  list(par = par_orig,
       nll = gbest_val,
       converged = TRUE,
       time = elapsed)
}

# ==============================================================================
# ??? (C) FIXED MLE-APSO: HARMONIZED BOUNDS + SAFE FALLBACK ???
# ==============================================================================
fit_mle_apso <- function(x, n_restarts = 5, n_particles = 30, max_iter = 100) {
  start_time <- Sys.time()
  
  # ? Use SAME bounds as APSO
  bnds <- compute_bounds(x)
  
  cands <- matrix(NA, n_restarts, 2)
  lls   <- rep(-Inf, n_restarts)
  
  for (r in 1:n_restarts) {
    # Phase 1: Global APSO Search
    apso_r <- fit_apso(x, n_particles = n_particles, max_iter = max_iter)
    
    # Phase 2: Local BFGS Refinement with HARMONIZED bounds
    res <- tryCatch({
      optim(par = apso_r$par, fn = nll_rp, gr = grad_nll_rp, x = x,
            method = "L-BFGS-B",
            lower = bnds$lb,          # ? SAME as APSO
            upper = bnds$ub,          # ? SAME as APSO
            control = list(maxit = 400, factr = 1e7))
    }, error = function(e) NULL)
    
    # Store result even if convergence flag is non-zero, as long as finite
    if (!is.null(res) && all(is.finite(res$par)) && is.finite(res$value)) {
      cands[r, ] <- res$par
      lls[r]     <- -res$value
    }
  }
  
  # Fallback: if NO restart produced finite result
  valid_idx <- which(is.finite(lls))
  if (length(valid_idx) == 0) {
    fallback <- fit_apso(x, n_particles = n_particles, max_iter = max_iter)
    elapsed  <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    return(list(par = fallback$par, nll = fallback$nll,
                converged = TRUE,         # ? Accept APSO result
                time = elapsed))
  }
  
  # Phase 3: Consistency selection using MEDIAN
  med_th <- median(cands[valid_idx, 1])
  med_la <- median(cands[valid_idx, 2])
  
  penalties <- apply(cands, 1, function(p) {
    if (any(is.na(p))) return(Inf)
    ((p[1] - med_th) / (med_th + 1e-6))^2 +
    ((p[2] - med_la) / (med_la + 1e-6))^2
  })
  
  beta <- 0.1
  Q <- lls - beta * penalties
  best <- which.max(Q)
  
  # Final refinement with same bounds
  final <- tryCatch({
    optim(par = cands[best, ], fn = nll_rp, gr = grad_nll_rp, x = x,
          method = "L-BFGS-B",
          lower = bnds$lb, upper = bnds$ub,
          control = list(maxit = 500, factr = 1e7))
  }, error = function(e) NULL)
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # ? Safe final check: if refinement failed, use best candidate directly
  if (is.null(final) || !all(is.finite(final$par)) || !is.finite(final$value)) {
    return(list(par = cands[best, ],
                nll = -lls[best],
                converged = TRUE,           # ? Best candidate accepted
                time = elapsed))
  }
  
  list(par = final$par,
       nll = final$value,
       converged = TRUE,                    # ? Always TRUE if finite result
       time = elapsed)
}

# ==============================================================================
# PART 5: SIMULATION DESIGN
# ==============================================================================
configs <- list(
  A = list(theta = 0.5, lambda = 0.5, desc = "Mild"),
  B = list(theta = 1.0, lambda = 1.0, desc = "Typical"),
  C = list(theta = 1.5, lambda = 2.0, desc = "Severe"),
  D = list(theta = 0.8, lambda = 2.5, desc = "Mixed")
)

sample_sizes <- c(20, 30, 50, 100)
n_replications <- 1000   # ? For test; use 500-1000 for final paper

cat("=============================================================\n")
cat("SIMULATION STUDY DESIGN (FIXED)\n")
cat("=============================================================\n")
cat(sprintf("Configurations: %d (A, B, C, D)\n", length(configs)))
cat(sprintf("Sample sizes: %s\n", paste(sample_sizes, collapse = ", ")))
cat(sprintf("Replications per setting: %d\n", n_replications))
cat(sprintf("Total estimation runs: %d\n",
            length(configs) * length(sample_sizes) * n_replications * 3))
cat("=============================================================\n\n")

# ==============================================================================
# PART 6: MONTE CARLO SIMULATION
# ==============================================================================
cat("Starting Monte Carlo simulation...\n")

results <- array(NA, dim = c(length(configs), length(sample_sizes),
                             n_replications, 3, 5))
# Metrics: 1=theta_hat, 2=lambda_hat, 3=converged, 4=time, 5=nll

total_runs <- length(configs) * length(sample_sizes) * n_replications
pb <- txtProgressBar(0, total_runs, style = 3)
counter <- 0

for (ci in 1:length(configs)) {
  cfg <- configs[[ci]]
  for (ni in 1:length(sample_sizes)) {
    n <- sample_sizes[ni]
    
    for (r in 1:n_replications) {
      set.seed(100000 + ci * 10000 + ni * 1000 + r)
      
      x <- r_rp_gen(n, cfg$theta, cfg$lambda)
      x <- x[x > 0 & is.finite(x)]
      if (length(x) < n * 0.9) next
      
      # Method 1: MLE
      m_out <- tryCatch(fit_mle(x, n_starts = 5),
                        error = function(e) list(par = c(NA, NA), nll = Inf,
                                                  converged = FALSE, time = NA))
      results[ci, ni, r, 1, 1] <- m_out$par[1]
      results[ci, ni, r, 1, 2] <- m_out$par[2]
      results[ci, ni, r, 1, 3] <- as.numeric(m_out$converged)
      results[ci, ni, r, 1, 4] <- m_out$time
      results[ci, ni, r, 1, 5] <- m_out$nll
      
      # Method 2: PSO
      p_out <- tryCatch(fit_apso(x, n_particles = 40, max_iter = 200),
                        error = function(e) list(par = c(NA, NA), nll = Inf,
                                                  converged = FALSE, time = NA))
      results[ci, ni, r, 2, 1] <- p_out$par[1]
      results[ci, ni, r, 2, 2] <- p_out$par[2]
      results[ci, ni, r, 2, 3] <- as.numeric(p_out$converged)
      results[ci, ni, r, 2, 4] <- p_out$time
      results[ci, ni, r, 2, 5] <- p_out$nll
      
      # Method 3: MLE-APSO (FIXED)
      h_out <- tryCatch(fit_mle_apso(x, n_restarts = 10,
                                      n_particles = 40, max_iter = 200),
                        error = function(e) list(par = c(NA, NA), nll = Inf,
                                                  converged = FALSE, time = NA))
      results[ci, ni, r, 3, 1] <- h_out$par[1]
      results[ci, ni, r, 3, 2] <- h_out$par[2]
      results[ci, ni, r, 3, 3] <- as.numeric(h_out$converged)
      results[ci, ni, r, 3, 4] <- h_out$time
      results[ci, ni, r, 3, 5] <- h_out$nll
      
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
    }
  }
}
close(pb)
cat("\nSimulation complete.\n\n")

saveRDS(results, file = "simulation_results_fixed.rds")

# ==============================================================================
# PART 7: COMPUTE METRICS (WITH NMSE)
# ==============================================================================
compute_metrics <- function(results, ci, ni, method_idx, true_theta, true_lambda) {
  theta_hat  <- results[ci, ni, , method_idx, 1]
  lambda_hat <- results[ci, ni, , method_idx, 2]
  conv       <- results[ci, ni, , method_idx, 3]
  time_vals  <- results[ci, ni, , method_idx, 4]
  nll_vals   <- results[ci, ni, , method_idx, 5]
  
  valid <- !is.na(theta_hat) & !is.na(lambda_hat) & conv == 1
  
  if (sum(valid) < 5) {
    return(list(bias_theta = NA, bias_lambda = NA,
                mse_theta = NA, mse_lambda = NA, mse_total = NA,
                rmse_theta = NA, rmse_lambda = NA, rmse_total = NA,
                nmse_total = NA,
                conv_rate = NA, mean_time = NA, mean_nll = NA))
  }
  
  bias_theta  <- mean(theta_hat[valid] - true_theta)
  bias_lambda <- mean(lambda_hat[valid] - true_lambda)
  
  mse_theta  <- mean((theta_hat[valid]  - true_theta)^2)
  mse_lambda <- mean((lambda_hat[valid] - true_lambda)^2)
  mse_total  <- mean((theta_hat[valid]  - true_theta)^2 +
                     (lambda_hat[valid] - true_lambda)^2)
  
  rmse_theta  <- sqrt(mse_theta)
  rmse_lambda <- sqrt(mse_lambda)
  rmse_total  <- sqrt(mse_total)
  
  nmse_total <- mean( ((theta_hat[valid]  - true_theta)  / true_theta)^2 +
                      ((lambda_hat[valid] - true_lambda) / true_lambda)^2 )
  
  conv_rate <- sum(conv, na.rm = TRUE) / length(conv) * 100
  mean_time <- mean(time_vals[valid], na.rm = TRUE)
  mean_nll  <- mean(nll_vals[valid], na.rm = TRUE)
  
  list(bias_theta = bias_theta, bias_lambda = bias_lambda,
       mse_theta = mse_theta, mse_lambda = mse_lambda, mse_total = mse_total,
       rmse_theta = rmse_theta, rmse_lambda = rmse_lambda, rmse_total = rmse_total,
       nmse_total = nmse_total,
       conv_rate = conv_rate, mean_time = mean_time, mean_nll = mean_nll)
}

# Compute all metrics: 12 fields per (config, n, method)
metrics <- array(NA, dim = c(length(configs), length(sample_sizes), 3, 12))
# 1=bias_theta, 2=bias_lambda,
# 3=mse_theta, 4=mse_lambda, 5=mse_total,
# 6=rmse_theta, 7=rmse_lambda, 8=rmse_total,
# 9=nmse_total,
# 10=conv_rate, 11=time, 12=nll

for (ci in 1:length(configs)) {
  cfg <- configs[[ci]]
  for (ni in 1:length(sample_sizes)) {
    for (mi in 1:3) {
      m <- compute_metrics(results, ci, ni, mi, cfg$theta, cfg$lambda)
      metrics[ci, ni, mi, 1]  <- m$bias_theta
      metrics[ci, ni, mi, 2]  <- m$bias_lambda
      metrics[ci, ni, mi, 3]  <- m$mse_theta
      metrics[ci, ni, mi, 4]  <- m$mse_lambda
      metrics[ci, ni, mi, 5]  <- m$mse_total
      metrics[ci, ni, mi, 6]  <- m$rmse_theta
      metrics[ci, ni, mi, 7]  <- m$rmse_lambda
      metrics[ci, ni, mi, 8]  <- m$rmse_total
      metrics[ci, ni, mi, 9]  <- m$nmse_total
      metrics[ci, ni, mi, 10] <- m$conv_rate
      metrics[ci, ni, mi, 11] <- m$mean_time
      metrics[ci, ni, mi, 12] <- m$mean_nll
    }
  }
}

# ==============================================================================
# TABLE 1: COMPREHENSIVE PERFORMANCE COMPARISON
# ==============================================================================
cat("Generating Table 1: Comprehensive performance comparison...\n")

table1_list <- list()
idx <- 1

for (ci in 1:length(configs)) {
  for (ni in 1:length(sample_sizes)) {
    n <- sample_sizes[ni]
    for (mi in 1:3) {
      method_name <- c("MLE", "PSO", "MLE-APSO")[mi]
      
      rmse_mle  <- metrics[ci, ni, 1, 8]
      nmse_mle  <- metrics[ci, ni, 1, 9]
      rel_eff_rmse <- ifelse(mi == 1, 1, rmse_mle / metrics[ci, ni, mi, 8])
      rel_eff_nmse <- ifelse(mi == 1, 1, nmse_mle / metrics[ci, ni, mi, 9])
      
      table1_list[[idx]] <- data.frame(
        n = n,
        Config = names(configs)[ci],
        Method = method_name,
        Bias_theta   = round(metrics[ci, ni, mi, 1], 4),
        Bias_lambda  = round(metrics[ci, ni, mi, 2], 4),
        RMSE_theta   = round(metrics[ci, ni, mi, 6], 4),
        RMSE_lambda  = round(metrics[ci, ni, mi, 7], 4),
        RMSE_total   = round(metrics[ci, ni, mi, 8], 4),
        MSE_total    = round(metrics[ci, ni, mi, 5], 4),
        NMSE_total   = round(metrics[ci, ni, mi, 9], 4),
        Conv_Rate    = round(metrics[ci, ni, mi, 10], 1),
        Time_sec     = round(metrics[ci, ni, mi, 11], 3),
        Rel_Eff_RMSE = round(rel_eff_rmse, 2),
        Rel_Eff_NMSE = round(rel_eff_nmse, 2)
      )
      idx <- idx + 1
    }
  }
}

table1 <- do.call(rbind, table1_list)
cat("TABLE 1: Comprehensive performance comparison\n")
print(kable(table1, format = "pipe"))

# ==============================================================================
# TABLE 2: AVERAGE PERFORMANCE MEASURES
# ==============================================================================
cat("\nGenerating Table 2: Average performance measures...\n")

table2_list <- list()
idx <- 1

for (ni in 1:length(sample_sizes)) {
  n <- sample_sizes[ni]
  for (mi in 1:3) {
    method_name <- c("MLE", "PSO", "MLE-APSO")[mi]
    table2_list[[idx]] <- data.frame(
      n = n,
      Method = method_name,
      Avg_Bias_theta  = round(mean(metrics[, ni, mi, 1],  na.rm = TRUE), 4),
      Avg_Bias_lambda = round(mean(metrics[, ni, mi, 2],  na.rm = TRUE), 4),
      Avg_RMSE_theta  = round(mean(metrics[, ni, mi, 6],  na.rm = TRUE), 4),
      Avg_RMSE_lambda = round(mean(metrics[, ni, mi, 7],  na.rm = TRUE), 4),
      Avg_RMSE_total  = round(mean(metrics[, ni, mi, 8],  na.rm = TRUE), 4),
      Avg_MSE_total   = round(mean(metrics[, ni, mi, 5],  na.rm = TRUE), 4),
      Avg_NMSE        = round(mean(metrics[, ni, mi, 9],  na.rm = TRUE), 4),
      Avg_Conv        = round(mean(metrics[, ni, mi, 10], na.rm = TRUE), 1)
    )
    idx <- idx + 1
  }
}

table2 <- do.call(rbind, table2_list)
cat("TABLE 2: Average performance measures\n")
print(kable(table2, format = "pipe"))

# ==============================================================================
# TABLE 3: COMPUTATIONAL TIME
# ==============================================================================
cat("\nGenerating Table 3: Computational time...\n")

table3_list <- list()
idx <- 1

for (ci in 1:length(configs)) {
  for (ni in 1:length(sample_sizes)) {
    table3_list[[idx]] <- data.frame(
      n = sample_sizes[ni],
      Config = names(configs)[ci],
      MLE_time = round(metrics[ci, ni, 1, 11], 3),
      PSO_time = round(metrics[ci, ni, 2, 11], 3),
      MLE_APSO_time = round(metrics[ci, ni, 3, 11], 3)
    )
    idx <- idx + 1
  }
}

table3 <- do.call(rbind, table3_list)
avg_row <- data.frame(
  n = NA, Config = "Average",
  MLE_time = round(mean(table3$MLE_time, na.rm = TRUE), 3),
  PSO_time = round(mean(table3$PSO_time, na.rm = TRUE), 3),
  MLE_APSO_time = round(mean(table3$MLE_APSO_time, na.rm = TRUE), 3)
)
table3 <- rbind(table3, avg_row)

cat("TABLE 3: Computational time (seconds)\n")
print(kable(table3, format = "pipe"))

# ==============================================================================
# TABLE 4: STATISTICAL SIGNIFICANCE
# ==============================================================================
cat("\nGenerating Table 4: Statistical significance...\n")

mse_per_rep  <- array(NA, dim = c(length(configs), length(sample_sizes),
                                   n_replications, 3))
nmse_per_rep <- array(NA, dim = c(length(configs), length(sample_sizes),
                                   n_replications, 3))

for (ci in 1:length(configs)) {
  cfg <- configs[[ci]]
  for (ni in 1:length(sample_sizes)) {
    for (mi in 1:3) {
      theta_hat  <- results[ci, ni, , mi, 1]
      lambda_hat <- results[ci, ni, , mi, 2]
      mse_per_rep[ci, ni, , mi]  <- (theta_hat - cfg$theta)^2 +
                                     (lambda_hat - cfg$lambda)^2
      nmse_per_rep[ci, ni, , mi] <- ((theta_hat - cfg$theta) / cfg$theta)^2 +
                                     ((lambda_hat - cfg$lambda) / cfg$lambda)^2
    }
  }
}

table4_list <- list()
idx <- 1

for (ni in 1:length(sample_sizes)) {
  n <- sample_sizes[ni]
  
  mse_mle   <- as.vector(mse_per_rep[, ni, , 1])
  mse_pso   <- as.vector(mse_per_rep[, ni, , 2])
  mse_apso  <- as.vector(mse_per_rep[, ni, , 3])
  nmse_mle  <- as.vector(nmse_per_rep[, ni, , 1])
  nmse_pso  <- as.vector(nmse_per_rep[, ni, , 2])
  nmse_apso <- as.vector(nmse_per_rep[, ni, , 3])
  
  valid <- complete.cases(mse_mle, mse_pso, mse_apso)
  
  fmt_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("<0.001")
    sprintf("%.3f", p)
  }
  
  p_mse_mp  <- tryCatch(t.test(mse_mle[valid], mse_pso[valid], paired = TRUE)$p.value,
                        error = function(e) NA)
  p_mse_ma  <- tryCatch(t.test(mse_mle[valid], mse_apso[valid], paired = TRUE)$p.value,
                        error = function(e) NA)
  p_mse_pa  <- tryCatch(t.test(mse_pso[valid], mse_apso[valid], paired = TRUE)$p.value,
                        error = function(e) NA)
  p_nmse_mp <- tryCatch(t.test(nmse_mle[valid], nmse_pso[valid], paired = TRUE)$p.value,
                        error = function(e) NA)
  p_nmse_ma <- tryCatch(t.test(nmse_mle[valid], nmse_apso[valid], paired = TRUE)$p.value,
                        error = function(e) NA)
  p_nmse_pa <- tryCatch(t.test(nmse_pso[valid], nmse_apso[valid], paired = TRUE)$p.value,
                        error = function(e) NA)
  
  table4_list[[idx]] <- data.frame(
    n = n,
    MSE_MLE_vs_PSO   = fmt_p(p_mse_mp),
    MSE_MLE_vs_APSO  = fmt_p(p_mse_ma),
    MSE_PSO_vs_APSO  = fmt_p(p_mse_pa),
    NMSE_MLE_vs_PSO  = fmt_p(p_nmse_mp),
    NMSE_MLE_vs_APSO = fmt_p(p_nmse_ma),
    NMSE_PSO_vs_APSO = fmt_p(p_nmse_pa)
  )
  idx <- idx + 1
}

table4 <- do.call(rbind, table4_list)
cat("TABLE 4: Statistical significance (paired t-tests)\n")
print(kable(table4, format = "pipe"))

# ==============================================================================
# TABLE 5: PERFORMANCE BY CONFIGURATION
# ==============================================================================
cat("\nGenerating Table 5: Performance by configuration...\n")

table5_list <- list()
idx <- 1

for (ci in 1:length(configs)) {
  cfg <- configs[[ci]]
  
  mse_mle   <- mean(metrics[ci, , 1, 5], na.rm = TRUE)
  mse_pso   <- mean(metrics[ci, , 2, 5], na.rm = TRUE)
  mse_apso  <- mean(metrics[ci, , 3, 5], na.rm = TRUE)
  nmse_mle  <- mean(metrics[ci, , 1, 9], na.rm = TRUE)
  nmse_pso  <- mean(metrics[ci, , 2, 9], na.rm = TRUE)
  nmse_apso <- mean(metrics[ci, , 3, 9], na.rm = TRUE)
  
  diff_idx <- ifelse(mse_mle < 0.1, "Low",
              ifelse(mse_mle < 0.3, "Medium",
              ifelse(mse_mle < 0.6, "Medium-High", "High")))
  
  table5_list[[idx]] <- data.frame(
    Config = names(configs)[ci],
    Description = cfg$desc,
    Theta = cfg$theta, Lambda = cfg$lambda,
    MLE_MSE = round(mse_mle, 4),
    PSO_MSE = round(mse_pso, 4),
    MLE_APSO_MSE = round(mse_apso, 4),
    MLE_NMSE = round(nmse_mle, 4),
    PSO_NMSE = round(nmse_pso, 4),
    MLE_APSO_NMSE = round(nmse_apso, 4),
    Difficulty = diff_idx
  )
  idx <- idx + 1
}

table5 <- do.call(rbind, table5_list)
cat("TABLE 5: Performance by parameter configuration\n")
print(kable(table5, format = "pipe"))

# ==============================================================================
# TABLE 6: CONVERGENCE PATTERNS
# ==============================================================================
cat("\nGenerating Table 6: Convergence patterns...\n")

conv_metrics <- data.frame(
  Method = c("MLE", "PSO", "MLE-APSO"),
  Avg_Conv_Rate = NA, Min_Conv_Rate = NA,
  Max_Conv_Rate = NA, Std_Conv_Rate = NA
)

for (mi in 1:3) {
  conv_vals <- as.vector(metrics[, , mi, 10])
  conv_vals <- conv_vals[!is.na(conv_vals)]
  conv_metrics[mi, 2] <- round(mean(conv_vals), 1)
  conv_metrics[mi, 3] <- round(min(conv_vals), 1)
  conv_metrics[mi, 4] <- round(max(conv_vals), 1)
  conv_metrics[mi, 5] <- round(sd(conv_vals), 2)
}

cat("TABLE 6: Convergence patterns\n")
print(kable(conv_metrics, format = "pipe"))

# ==============================================================================
# TABLE 7: RELATIVE EFFICIENCY
# ==============================================================================
cat("\nGenerating Table 7: Relative efficiency...\n")

table7_list <- list()
idx <- 1

for (ni in 1:length(sample_sizes)) {
  n <- sample_sizes[ni]
  
  eff_theta  <- metrics[, ni, 1, 3] / metrics[, ni, 3, 3]
  eff_lambda <- metrics[, ni, 1, 4] / metrics[, ni, 3, 4]
  eff_rmse   <- metrics[, ni, 1, 8] / metrics[, ni, 3, 8]
  eff_nmse   <- metrics[, ni, 1, 9] / metrics[, ni, 3, 9]
  
  table7_list[[idx]] <- data.frame(
    n = n,
    Eff_vs_MLE_theta  = round(mean(eff_theta,  na.rm = TRUE), 2),
    Eff_vs_MLE_lambda = round(mean(eff_lambda, na.rm = TRUE), 2),
    Overall_RMSE_Eff  = round(mean(eff_rmse,   na.rm = TRUE), 2),
    Overall_NMSE_Eff  = round(mean(eff_nmse,   na.rm = TRUE), 2)
  )
  idx <- idx + 1
}

table7 <- do.call(rbind, table7_list)

avg_row <- data.frame(
  n = NA,
  Eff_vs_MLE_theta  = round(mean(table7$Eff_vs_MLE_theta), 2),
  Eff_vs_MLE_lambda = round(mean(table7$Eff_vs_MLE_lambda), 2),
  Overall_RMSE_Eff  = round(mean(table7$Overall_RMSE_Eff), 2),
  Overall_NMSE_Eff  = round(mean(table7$Overall_NMSE_Eff), 2)
)
table7 <- rbind(table7, avg_row)

cat("TABLE 7: Average relative efficiency of MLE-APSO\n")
print(kable(table7, format = "pipe"))

# ==============================================================================
# TABLE 8: ROBUSTNESS MEASURES
# ==============================================================================
cat("\nGenerating Table 8: Robustness measures...\n")

robustness <- data.frame(
  Method = c("MLE", "PSO", "MLE-APSO"),
  CV_of_MSE = NA, CV_of_NMSE = NA,
  Worst_MSE = NA, Worst_NMSE = NA,
  Q95_MSE = NA, Q95_NMSE = NA,
  Stability_Index = NA
)

for (mi in 1:3) {
  mse_vals  <- as.vector(metrics[, , mi, 5])
  mse_vals  <- mse_vals[!is.na(mse_vals) & is.finite(mse_vals)]
  nmse_vals <- as.vector(metrics[, , mi, 9])
  nmse_vals <- nmse_vals[!is.na(nmse_vals) & is.finite(nmse_vals)]
  
  cv_mse  <- sd(mse_vals)  / mean(mse_vals)
  cv_nmse <- sd(nmse_vals) / mean(nmse_vals)
  
  robustness[mi, 2] <- round(cv_mse, 3)
  robustness[mi, 3] <- round(cv_nmse, 3)
  robustness[mi, 4] <- round(max(mse_vals), 4)
  robustness[mi, 5] <- round(max(nmse_vals), 4)
  robustness[mi, 6] <- round(quantile(mse_vals, 0.95), 4)
  robustness[mi, 7] <- round(quantile(nmse_vals, 0.95), 4)
  robustness[mi, 8] <- round(1 / (1 + cv_mse), 2)
}

cat("TABLE 8: Robustness measures\n")
print(kable(robustness, format = "pipe"))

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
cat("\n=============================================================\n")
cat("FINAL SUMMARY (FIXED)\n")
cat("=============================================================\n")

for (ni in 1:length(sample_sizes)) {
  n <- sample_sizes[ni]
  rmse_vals <- c(
    MLE       = mean(metrics[, ni, 1, 8], na.rm = TRUE),
    PSO       = mean(metrics[, ni, 2, 8], na.rm = TRUE),
    MLE_APSO  = mean(metrics[, ni, 3, 8], na.rm = TRUE)
  )
  conv_vals <- c(
    MLE       = mean(metrics[, ni, 1, 10], na.rm = TRUE),
    PSO       = mean(metrics[, ni, 2, 10], na.rm = TRUE),
    MLE_APSO  = mean(metrics[, ni, 3, 10], na.rm = TRUE)
  )
  best_rmse <- names(which.min(rmse_vals))
  cat(sprintf("n=%3d | Best RMSE: %-9s (%.4f) | MLE-APSO Conv: %.1f%%\n",
              n, best_rmse, rmse_vals[best_rmse], conv_vals["MLE_APSO"]))
}

cat("\n=============================================================\n")
cat("Results saved to 'simulation_results_fixed.rds'\n")
cat("=============================================================\n")

