## R package 설치
packages <- c("tidyverse", "rstan", "parallel")

installed <- rownames(installed.packages())
for(p in packages){if(!(p %in% installed)){install.packages(p)}}

invisible(lapply(packages, library, character.only = TRUE))

source("LTBI_data_description.r")

## 데이터 관련 설정
data_start <- 2017
data_end <- 2024
age <- 19

years_all <- (data_start - age):data_end
T_all <- length(years_all)


## 모델 파라미터 관련 설정
## reversion
rev_annual_prop <- 0.01
gamma_est <- -log(1 - rev_annual_prop)

## 출생 시 감수성 분율
prev_0 <- 1


## 추가 데이터 활용: 고1 대상 LTBI 양성률 데이터 (IGRA 기반)
year <- c(2015, 2017, 2018)
age_obs <- 15
sample <- c(58431, 270435, 20587)
positive <- c(1029, 5476, 373)

df_obs <- bind_rows(
    df_est %>% transmute(year = as.integer(year), age_obs = 19, 
                         sample = as.integer(sample), positive = as.integer(positive), prev = positive / sample),
    tibble(year, age_obs, sample, positive) %>% mutate(prev = positive / sample)) %>% arrange(year, age_obs)

T_obs <- nrow(df_obs)
t_obs_idx <- match(df_obs$year, years_all)
N_vec <- df_obs$sample
y_vec <- df_obs$positive


## IGRA 검사 관련
sens_mean_row <- pick_by_year(df_obs$year, sens_4th, sens_SD, sens_3rd)
spec_mean_row <- pick_by_year(df_obs$year, spec_4th, spec_SD, spec_3rd)

sens_sd_row <- pick_by_year(df_obs$year, logit_sd_from_ci(sens_ci_4th), logit_sd_from_ci(sens_ci_SD), logit_sd_from_ci(sens_ci_3rd))
spec_sd_row <- pick_by_year(df_obs$year, logit_sd_from_ci(spec_ci_4th), logit_sd_from_ci(spec_ci_SD), logit_sd_from_ci(spec_ci_3rd))

is15 <- df_obs$age_obs == 15
sens_mean_row[is15] <- sens_3rd
spec_mean_row[is15] <- spec_3rd
sens_sd_row[is15] <- logit_sd_from_ci(sens_ci_3rd)
spec_sd_row[is15] <- logit_sd_from_ci(spec_ci_3rd)

sens_mu_logit_row <- qlogis(pmin(pmax(sens_mean_row, 1e-9), 1 - 1e-9))
spec_mu_logit_row <- qlogis(pmin(pmax(spec_mean_row, 1e-9), 1 - 1e-9))
sens_sd_logit_row <- pmax(sens_sd_row, 1e-4)
spec_sd_logit_row <- pmax(spec_sd_row, 1e-4)


set.seed(1)
n_sim <- 10000
M_use <- 10

cl <- makeCluster(max(1L, detectCores() - 1L))
clusterExport(cl, c("df_obs", "sens_mu_logit_row", "sens_sd_logit_row", "spec_mu_logit_row",
                    "spec_sd_logit_row", "n_sim", "M_use", "rlogitnorm", "sim_one_row_RG_fast"), envir = environment())
clusterSetRNGStream(cl, 1)

res_list <- parLapplyLB( cl, seq_len(T_obs),
                        function(i) sim_one_row_RG_fast(
                            i, df_obs, sens_mu_logit_row, sens_sd_logit_row, spec_mu_logit_row, spec_sd_logit_row, n_sim, M_use))
stopCluster(cl)

sens_logit_feas <- qlogis(do.call(rbind, lapply(res_list, `[[`, "se")))
spec_logit_feas <- qlogis(do.call(rbind, lapply(res_list, `[[`, "sp")))


## Anchor 데이터 활용
## Lee et al., 2021, BMC ID
temp_ref <- tibble(year = c(1999, 2001, 2004, 2005),
                   med_ARI = c(.51, .36, .26, .29), 
                   upper_ARI = c(.62, .47, .34, .44), 
                   lower_ARI = c(.40, .26, .18, .15))

## reversion 관련 보정
adj_val <- 1.5

## 1999년 추정치를 1998년으로 가정
temp_ref2 <- bind_rows(tibble(year = 1998, 
                              med_ARI = temp_ref$med_ARI[temp_ref$year == 1999],
                              lower_ARI = temp_ref$lower_ARI[temp_ref$year == 1999],
                              upper_ARI = temp_ref$upper_ARI[temp_ref$year == 1999]), temp_ref)

ari_med <- pmin(pmax(temp_ref2$med_ARI * adj_val / 100, 1e-9), 1 - 1e-9)
ari_lo  <- pmin(pmax(temp_ref2$lower_ARI * adj_val / 100, 1e-9), 1 - 1e-9)
ari_hi  <- pmin(pmax(temp_ref2$upper_ARI * adj_val / 100, 1e-9), 1 - 1e-9)

## ARI에서 force of infection으로 변환
lambda_med <- -log(1 - ari_med)
lambda_lo <- -log(1 - ari_lo)
lambda_hi <- -log(1 - ari_hi)

anchor_loglam_mu <- log(pmax(lambda_med, 1e-12))
anchor_loglam_sd <- pmax((log(pmax(lambda_hi, 1e-12)) - log(pmax(lambda_lo, 1e-12))) / (2 * 1.96), 0.5)
anchor_loglam_sd[temp_ref2$year == 1998] <- 1.3
anchor_idx <- match(temp_ref2$year, years_all)


## 연도별 ARI 추정: R stan을 활용한 베이지언 추정법
options(mc.cores = detectCores())
rstan_options(auto_write = TRUE)

adjusted_ARI_model <- "
data {
  int<lower=1> T;
  int<lower=1> T_obs;
  array[T_obs] int<lower=0> N;
  array[T_obs] int<lower=0> y;
  array[T_obs] int<lower=1, upper=T> t_obs_idx;
  array[T_obs] int<lower=1> age_obs;
  array[T_obs] int<lower=0, upper=1> use_row;
  real<lower=0> gamma;
  real<lower=0, upper=1> prev0;
  int<lower=1> M;
  matrix[T_obs, M] sens_logit_feas;
  matrix[T_obs, M] spec_logit_feas;
  real<lower=0> lambda0;
  int<lower=0> K_anchor;
  array[K_anchor] int<lower=1, upper=T> anchor_idx;
  vector[K_anchor] anchor_loglam_mu;
  vector<lower=0>[K_anchor] anchor_loglam_sd;
}
parameters {
  vector[T] log_lambda;
  real<lower=0> sigma_rw2;
  real<lower=0> phi;
}
transformed parameters {
  vector[T] lambda = exp(log_lambda);
  vector[T] p_true15;
  vector[T] p_true19;

  for (t in 1:T) {
    {
      int birth_idx = t - 19;
      real S = prev0;
      real I = 0;
      if (birth_idx < 1) p_true19[t] = 0;
      else {
        for (k in birth_idx:t) {
          real p_inf = 1 - exp(-lambda[k]);
          real newInf = S * p_inf;
          I = I * exp(-gamma) + newInf;
          S = S * exp(-lambda[k]);
        }
        p_true19[t] = fmin(fmax(I, 1e-12), 1 - 1e-12);
      }
    }
    {
      int birth_idx = t - 15;
      real S = prev0;
      real I = 0;
      if (birth_idx < 1) p_true15[t] = 0;
      else {
        for (k in birth_idx:t) {
          real p_inf = 1 - exp(-lambda[k]);
          real newInf = S * p_inf;
          I = I * exp(-gamma) + newInf;
          S = S * exp(-lambda[k]);
        }
        p_true15[t] = fmin(fmax(I, 1e-12), 1 - 1e-12);
      }
    }
  }
}
model {
  log_lambda[1] ~ normal(log(lambda0), 1);
  log_lambda[2] ~ normal(log_lambda[1], 1);
  for (t in 3:T) log_lambda[t] ~ normal(2 * log_lambda[t - 1] - log_lambda[t - 2], sigma_rw2);
  sigma_rw2 ~ normal(0, 0.2) T[0,];
  phi ~ lognormal(log(200), 0.7);

  for (k in 1:K_anchor)
    log_lambda[anchor_idx[k]] ~ normal(anchor_loglam_mu[k], anchor_loglam_sd[k]);

  for (n in 1:T_obs) if (use_row[n] == 1) {
    int tt = t_obs_idx[n];
    real p_true_use = age_obs[n] == 19 ? p_true19[tt] : p_true15[tt];
    vector[M] se = inv_logit(to_vector(sens_logit_feas[n]));
    vector[M] sp = inv_logit(to_vector(spec_logit_feas[n]));
    vector[M] lp;

    for (m in 1:M) {
      real p_obs = se[m] * p_true_use + (1 - sp[m]) * (1 - p_true_use);
      real a = fmin(fmax(p_obs, 1e-12), 1 - 1e-12) * phi;
      real b = (1 - fmin(fmax(p_obs, 1e-12), 1 - 1e-12)) * phi;
      lp[m] = beta_binomial_lpmf(y[n] | N[n], a, b);
    }
    target += log_sum_exp(lp) - log(M);
  }
}
generated quantities {
  vector[T] ARI;
  for (t in 1:T) ARI[t] = 1 - exp(-exp(log_lambda[t]));
}
"


use_row <- rep(1L, T_obs)
use_row[df_obs$age_obs == 19L & df_obs$year %in% c(2021L, 2023L)] <- 0L
lambda0_val <- 0.008

stan_data <- list(T = T_all, T_obs = T_obs, N = as.integer(N_vec), y = as.integer(y_vec), t_obs_idx = as.integer(t_obs_idx),
                  age_obs = as.integer(df_obs$age_obs), use_row = as.integer(use_row),
                  gamma = gamma_est, prev0 = prev_0, M = M_use,
                  sens_logit_feas = sens_logit_feas, spec_logit_feas = spec_logit_feas, lambda0 = lambda0_val, 
                  K_anchor = length(anchor_idx), anchor_idx = as.integer(anchor_idx), 
                  anchor_loglam_mu = as.numeric(anchor_loglam_mu), anchor_loglam_sd = as.numeric(anchor_loglam_sd))

fit <- stan(model_code = adjusted_ARI_model, data = stan_data, chains = 4, iter = 1000, warmup = 200,
            control = list(adapt_delta = 0.99, max_treedepth = 18),
            init = function() list(log_lambda = rnorm(T_all, log(lambda0_val), 0.6),
                                   sigma_rw2 = abs(rnorm(1, 0.04, 0.02))))


## 추정된 ARI 시각화
ARI_draws <- as.matrix(extract(fit, "ARI")$ARI)

est_ARI <- tibble(year = years_all, med_ARI = apply(ARI_draws, 2, median) * 100,
                  lower_ARI = apply(ARI_draws, 2, quantile, 0.025) * 100,
                  upper_ARI = apply(ARI_draws, 2, quantile, 0.975) * 100)

ggplot(est_ARI, aes(x = year, y = med_ARI)) +
geom_line(linewidth = 1.2) +
geom_ribbon(aes(ymin = lower_ARI, ymax = upper_ARI), alpha = 0.25) +
geom_point(data = temp_ref, aes(x = year, y = med_ARI*adj_val), color = "blue", size = 3) +
geom_errorbar(data = temp_ref, aes(x = year, ymin = lower_ARI*adj_val, ymax = upper_ARI*adj_val), width = 0.3, color = "blue") +
labs(x = "Year", y = "Time-varying annual risk of LTBI (%)") +
scale_x_continuous(expand = c(0, 0), breaks = c(2000, 2005, 2010, 2015, 2020)) +
scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0)) +
theme(axis.title = element_text(size = 17, family = "sans", colour = "black"),
      axis.text = element_text(size = 15, family = "sans", color = "black"))



## LTBI 유병률 추정을 위한 과거로의 ARI 외삽
## 외삽 시작 연도 설정
year_extra <- 1920

## TST 기반 ARI 데이터
tst_year <- c(1970, 1975, 1980, 1985, 1990)
tst_ari  <- c(3.9, 2.3, 1.8, 1.2, 1.1) / 100
tst_df <- data.frame(tst_year, tst_ari)

eps <- 1e-8
set.seed(1)

years_fit <- years_all
y0 <- min(years_fit)
years_pre <- year_extra:(y0 - 1L)
years_ext <- year_extra:max(years_fit)

ARI_draws <- as.matrix(rstan::extract(fit, pars = "ARI")$ARI)
keep_id <- sample(seq_len(nrow(ARI_draws)), min(2000, nrow(ARI_draws)))

x_fit <- years_fit - y0
u_pre <- abs(years_pre - y0)
tst_u <- abs(tst_year - y0)
tst_y <- qlogis(pmin(pmax(tst_ari, eps), 1 - eps))

y0_draw <- qlogis(pmin(pmax(ARI_draws[keep_id, 1], eps), 1 - eps))
y0_med <- median(y0_draw)

local_slope_post <- function(p_fit, win = 7, eps = 1e-8) {
    y <- qlogis(pmin(pmax(as.numeric(p_fit), eps), 1 - eps))
    k <- min(win, length(y) - 1L)
    s <- coef(lm(y[1:(k + 1)] ~ x_fit[1:(k + 1)]))[2]
    ifelse(is.finite(s), s, 0)
}

h_u <- function(u, tau) (1 - exp(-u / tau))^2
dh_u <- function(u, tau) 2 * (1 - exp(-u / tau)) * exp(-u / tau) / tau

win_slope <- 7
sd_tst <- 1
tst_w <- c(0.10, 1, 1, 1, 1)
rho <- 0.45
tau <- 25
sd_c1 <- 0.02
s_ref <- 0.01
sd_logs <- 0.35
cabs_ref <- 5e-4
sd_logcabs <- 0.60
mono_sd <- 0.005

s_post_typ <- median(vapply(keep_id, function(i) local_slope_post(ARI_draws[i, ], win_slope, eps), 0), na.rm = TRUE)

obj_sc <- function(theta) {
    s <- exp(theta[1])
    c <- -exp(theta[2])
    dev0 <- y0_draw - y0_med
    
    y_pred <- outer(y0_draw, s * tst_u + c * tst_u^2, "+") + (-rho) * outer(dev0, h_u(tst_u, tau), "*")
    res <- sweep(y_pred, 2, tst_y, "-")
    
    deriv <- s + 2 * c * u_pre + (-rho) * quantile(abs(dev0), 0.975, na.rm = TRUE) * dh_u(u_pre, tau)
    
    sum(tst_w * colSums((res / sd_tst)^2)) + ((-s - s_post_typ) / sd_c1)^2 + sum((pmax(0, -deriv) / mono_sd)^2) +
    ((theta[1] - log(s_ref)) / sd_logs)^2 + ((theta[2] - log(cabs_ref)) / sd_logcabs)^2
}

opt <- optim(c(log(max(1e-6, -s_post_typ)), log(cabs_ref)), obj_sc, method = "L-BFGS-B",
             lower = c(log(1e-6), log(1e-8)), upper = c(log(0.04), log(5e-2)),control = list(maxit = 800))

s_global <- exp(opt$par[1])
c_global <- -exp(opt$par[2])

extrap_with_sc <- function(p_fit, eps = 1e-8) {
    y0_1998 <- qlogis(pmin(pmax(as.numeric(p_fit[1]), eps), 1 - eps))
    delta_u <- (-rho) * (y0_1998 - y0_med) * h_u(u_pre, tau)
    c(plogis(y0_1998 + s_global * u_pre + c_global * u_pre^2 + delta_u), p_fit)
}

ari_ext_mat <- vapply(keep_id, function(i) extrap_with_sc(ARI_draws[i, ], eps), numeric(length(years_ext)))

ari_ext_df <- tibble(year = years_ext, 
                     med_ARI = apply(ari_ext_mat, 1, median) * 100,
                     lower_ARI = apply(ari_ext_mat, 1, quantile, 0.025) * 100,
                     upper_ARI = apply(ari_ext_mat, 1, quantile, 0.975) * 100)


## 외삽 ARI 시각화
ggplot(ari_ext_df, aes(year, med_ARI)) +
geom_ribbon(data = filter(ari_ext_df, year <= 1998), aes(ymin = lower_ARI, ymax = upper_ARI), fill = "blue", alpha = 0.25) +
geom_line(data = filter(ari_ext_df, year <= 1998), color = "blue", linewidth = 1) +
geom_ribbon(data = filter(ari_ext_df, year >= 1998), aes(ymin = lower_ARI, ymax = upper_ARI), alpha = 0.25) +
geom_line(data = filter(ari_ext_df, year >= 1998), linewidth = 1) +
geom_vline(xintercept = y0, linetype = 2) +
geom_point(data = tst_df, aes(tst_year, tst_ari*100), size = 2.5, inherit.aes = FALSE) +
labs(x = "Year", y = "Annual risk of LTBI (%)") +
scale_x_continuous(expand = c(0, 0)) +
scale_y_continuous(limits = c(0, 25), expand = c(0, 0)) +
theme(axis.title = element_text(size = 17, family = "sans", colour = "black"),
      axis.text = element_text(size = 15, family = "sans", color = "black"))


## 연령그룹별 LTBI 유병률 계산
## 인구 데이터 (0세부터 85+세까지 1세별 인구 데이터)
pop_age <- read.csv("../data/data_pop.csv")$pop

## 대상 연도
year_target <- 2024

years_ext <- as.integer(years_ext)
if (nrow(ari_ext_mat) != length(years_ext) && ncol(ari_ext_mat) == length(years_ext)) {ari_ext_mat <- t(ari_ext_mat)}
stopifnot(nrow(ari_ext_mat) == length(years_ext))

lambda_ext_mat <- -log(1 - pmin(pmax(ari_ext_mat, 1e-12), 1 - 1e-12))

gamma <- 0
t_eval <- 2024L
ages_target <- 0:85

pop_age <- read.csv("../data/data_pop.csv")$pop
pop_age[!is.finite(pop_age)] <- 0
if (sum(pop_age) <= 0) stop("pop_age sums to 0")

prev_from_lambda_path <- function(lambda_path_mat, gamma, prev0 = 1) {
    S <- rep(prev0, ncol(lambda_path_mat))
    I <- rep(0, ncol(lambda_path_mat))
    for (j in seq_len(nrow(lambda_path_mat))) {
        p_inf <- 1 - exp(-lambda_path_mat[j, ])
        newInf <- S * p_inf
        I <- I * exp(-gamma) + newInf
        S <- S * exp(-lambda_path_mat[j, ])
    }
    pmin(pmax(I, 1e-12), 1 - 1e-12)
}

prev_draws_2024 <- sapply(ages_target, function(a) {
    idx <- match((t_eval - a):t_eval, years_ext)
    if (any(is.na(idx))) stop("years_ext does not cover age ", a)
    prev_from_lambda_path(lambda_ext_mat[idx, , drop = FALSE], gamma)}) %>% t()

age_group <- cut(ages_target, breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, Inf),
    right = FALSE,
    labels = c("0–4세","5–9세","10–14세","15–19세","20–24세","25–29세","30–34세",
               "35–39세","40–44세","45–49세","50–54세","55–59세","60–64세","65세 이상"))

wmean_draws <- function(x, w) as.numeric(crossprod(w / sum(w), x))
summarise_draws <- function(x) c(med = median(x), lo = quantile(x, 0.025), hi = quantile(x, 0.975))
fmt_ci <- function(x) sprintf("%.2f%% (%.2f–%.2f%%)", 100 * x[1], 100 * x[2], 100 * x[3])

group_draws <- lapply(levels(age_group), function(g) {
    ii <- which(age_group == g)
    wmean_draws(prev_draws_2024[ii, , drop = FALSE], pop_age[ii])
})
names(group_draws) <- levels(age_group)
all_draws <- wmean_draws(prev_draws_2024, pop_age)

prev_table_2024 <- bind_rows(
    lapply(names(group_draws), function(g) {
        tibble(연령군 = g, `유병률(95% CI)` = fmt_ci(summarise_draws(group_draws[[g]])))}),
    tibble(연령군 = "전 연령", `유병률(95% CI)` = fmt_ci(summarise_draws(all_draws)))
)

print(prev_table_2024, n = Inf)