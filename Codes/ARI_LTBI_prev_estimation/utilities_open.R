
pick_by_year <- function(year_vec, v_4th, v_SD, v_3rd) {
    out <- rep(v_3rd, length(year_vec))
    out[year_vec %in% years_3rd] <- v_3rd
    out[year_vec %in% years_4th] <- v_4th
    out[year_vec %in% years_SD]  <- v_SD
    return(out)
}




logit_sd_from_ci <- function(ci, eps = 1e-6) {
    L <- pmin(pmax(ci[1], eps), 1 - eps)
    U <- pmin(pmax(ci[2], eps), 1 - eps)
    return((qlogis(U) - qlogis(L)) / (2*1.96))
}




sim_one_row_RG_new <- function(i, df, sens_mean, sens_sd, spec_mean, spec_sd, n_sim) { 
    mu_s <- qlogis(sens_mean[i]); sd_s <- sens_sd[i]
    mu_c <- qlogis(spec_mean[i]); sd_c <- spec_sd[i]    
    s_draw <- plogis(rnorm(n_sim, mu_s, sd_s))
    c_draw <- plogis(rnorm(n_sim, mu_c, sd_c))
    list(sens = s_draw, spec = c_draw)
}




binom_wilson_ci <- function(x, n, conf.level = 0.95) {
    p_hat <- x / n
    z <- qnorm(1 - (1 - conf.level) / 2)
    denom <- 1 + z^2 / n
    centre <- (p_hat + z^2 / (2 * n)) / denom
    margin <- z * sqrt((p_hat * (1 - p_hat) + z^2 / (4 * n)) / n) / denom
    lower <- centre - margin
    upper <- centre + margin
    return(c(estimate = p_hat, lower = lower, upper = upper))
}




rlogitnorm <- function(n, mu_logit, sd_logit) plogis(rnorm(n, mu_logit, sd_logit))




sim_one_row_RG_fast <- function(i, df_obs, mu_se, sd_se, mu_sp, sd_sp, n_sim, M_use) {
    se <- rlogitnorm(n_sim, mu_se[i], sd_se[i])
    sp <- rlogitnorm(n_sim, mu_sp[i], sd_sp[i])
    den_ok <- (se + sp - 1) > 0
    
    ysim <- rbinom(n_sim, size = df_obs$sample[i], prob = df_obs$prev[i])
    prev_sim <- ysim / df_obs$sample[i]
    
    feas <- den_ok & (prev_sim >= (1 - sp))
    rate <- mean(feas)
    
    idx <- which(feas)
    pick <- sample(idx, size = M_use, replace = TRUE)
    list(se = se[pick], sp = sp[pick], feas_rate = rate)
}