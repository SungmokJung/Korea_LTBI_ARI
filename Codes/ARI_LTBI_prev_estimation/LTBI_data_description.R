## R package 설치
packages <- c("dplyr", "tidyr", "ggplot2", "parallel", "scales")

installed <- rownames(installed.packages())
for(p in packages){if(!(p %in% installed)){install.packages(p)}}

invisible(lapply(packages, library, character.only = TRUE))


## 함수 불러오기
source("utilities_open.r")


## 관찰된 LTBI 유병률 데이터: 병무청 IGRA 자료
## 관찰된 LTBI 유병률 데이터
data.frame(year = 2017:2024,
           positive = c(9731, 8347, 4351, 3049, 5094, 2568, 3775, 2482),
           sample = c(333739, 323771, 332526, 288814, 260104, 255231, 244559, 227537)) %>% 
mutate(prev = positive / sample) -> observed_data


## IGRA 검사 관련 데이터
## 연도별 IGRA 검사 정보
years_3rd <- c(2017, 2018)
years_4th <- c(2019, 2020, 2022, 2024)
years_SD  <- c(2021, 2023)

## IGRA 검사별 민감도 & 특이도
sens_3rd <- 0.81; sens_ci_3rd <- c(0.79, 0.84)
spec_3rd <- 0.99; spec_ci_3rd <- c(0.98, 0.99)

sens_4th <- 0.89; sens_ci_4th <- c(0.84, 0.94)
spec_4th <- 0.99; spec_ci_4th <- c(0.98, 0.99)

sens_SD <- 0.88; sens_ci_SD <- c(0.73, 0.96)
spec_SD <- 0.95; spec_ci_SD <- c(0.91, 0.98)


## 실제 유병률 및 양성예측도 추정: Rogan–Gladen 추정법
n_sim <- 50000
n_row <- nrow(observed_data)

## 민감도와 특이도 샘플링
sens_mean <- pick_by_year(observed_data$year, sens_4th, sens_SD, sens_3rd)
sens_sd <- pick_by_year(observed_data$year, logit_sd_from_ci(sens_ci_4th), logit_sd_from_ci(sens_ci_SD), logit_sd_from_ci(sens_ci_3rd))

spec_mean <- pick_by_year(observed_data$year, spec_4th, spec_SD, spec_3rd)
spec_sd <- pick_by_year(observed_data$year, logit_sd_from_ci(spec_ci_4th), logit_sd_from_ci(spec_ci_SD), logit_sd_from_ci(spec_ci_3rd))

cores_to_use <- max(1L, detectCores() - 1L)

mclapply(X = seq_len(n_row), FUN = sim_one_row_RG_new, df = observed_data,
         sens_mean = sens_mean, sens_sd = sens_sd, spec_mean = spec_mean, spec_sd = spec_sd, 
         n_sim = n_sim, mc.cores = cores_to_use) -> res_list

sens_mat <- do.call(rbind, lapply(res_list, `[[`, "sens"))
spec_mat <- do.call(rbind, lapply(res_list, `[[`, "spec"))

## 실제 양성률 및 양성예측도 추정
y_sim <- matrix(rbinom(n_row * n_sim, size = rep(observed_data$sample, each = n_sim), prob = rep(observed_data$prev, each = n_sim)),
                nrow = n_row, ncol = n_sim, byrow = TRUE)
prev_mat <- y_sim / matrix(observed_data$sample, nrow = n_row, ncol = n_sim)

den_raw <- sens_mat + spec_mat - 1
feasible <- (den_raw > 0) & (prev_mat >= (1 - spec_mat))

p_true_mat <- matrix(NA_real_, nrow = n_row, ncol = n_sim)
p_true_mat[feasible] <- (prev_mat[feasible] - (1 - spec_mat[feasible])) / den_raw[feasible]
p_true_mat <- pmin(pmax(p_true_mat, 0), 1)

ppv_mat <- matrix(NA_real_, nrow = n_row, ncol = n_sim)
ppv_mat[feasible] <- (sens_mat[feasible] * p_true_mat[feasible]) / 
(sens_mat[feasible] * p_true_mat[feasible] + (1 - spec_mat[feasible]) * (1 - p_true_mat[feasible]))


## 추정된 실제 양성률 및 양성예측도
tibble(year = observed_data$year,
       prev_true_med = apply(p_true_mat, 1, quantile, .5, na.rm = TRUE),
       prev_true_lower = apply(p_true_mat, 1, quantile, .025, na.rm = TRUE),
       prev_true_upper = apply(p_true_mat, 1, quantile, .975, na.rm = TRUE)) -> prev_summary

tibble(year = observed_data$year, 
       ppv_med = apply(ppv_mat, 1, quantile, .5, na.rm = TRUE), 
       ppv_lower = apply(ppv_mat, 1, quantile, .025, na.rm = TRUE),
       ppv_upper = apply(ppv_mat, 1, quantile, .975, na.rm = TRUE)) -> ppv_summary

observed_data %>% dplyr::select(year, sample, positive, prev) %>% 
left_join(prev_summary, by = "year") %>% left_join(ppv_summary, by = "year") -> df_est



## Rogan-Gladen 추정법 기반 이론적으로 가능한 sample
feas_rate <- data.frame(year = observed_data$year, feasible_prop = rowMeans(feasible))


## LTBI 양성률 시각화 
df_est %>% dplyr::select(year, sample, positive, prev) %>% rowwise() %>%
mutate(lower = binom_wilson_ci(positive, sample)[2], upper = binom_wilson_ci(positive, sample)[3]) -> df_ci

df_est %>% dplyr::select(year, prev, prev_true_med) %>% 
tidyr::pivot_longer(cols = -year, names_to = "group", values_to = "value") %>%
mutate(group = dplyr::case_when(group == "prev" ~ "Observed", group == "prev_true_med" ~ "Estimated")) -> df_figure

df_figure$group <- factor(df_figure$group, levels = c("Observed", "Estimated"))

options(repr.plot.width = 8, repr.plot.height = 6, warn = -1)
offset <- .22
theme_set(theme_bw())

ggplot(data = df_figure, aes(x = factor(year))) +
geom_bar(aes(y = value, group = group, fill = group), stat = "identity", position = 'dodge') +
geom_errorbar(data = df_ci, aes(x = as.numeric(factor(year)) - offset, ymin = lower, ymax = upper), width = .3) +
geom_errorbar(data = df_est, aes(x = as.numeric(factor(year)) + offset, ymin = prev_true_lower, ymax = prev_true_upper), width = .3) +
scale_y_continuous(limits = c(0, .04), expand = c(0, 0), labels = percent_format(accuracy = 1)) +
labs(x = "Year", y = "IGRA positivity") +
theme(axis.title = element_text(size = 17, family = "sans", colour = "black"),
      axis.text = element_text(size = 15, family = "sans", color = "black"),
      legend.title = element_blank(), legend.text = element_text(size = 15, color = "black"),
      legend.position = "inside", legend.position.inside = c(.85, .91))

## 양성예측도 시각화 
df_ppv_long <- as.data.frame(t(ppv_mat)) %>% 
mutate(draw = row_number()) %>% 
pivot_longer(-draw, names_to = "row_id", values_to = "ppv") %>%
mutate(row_id = as.integer(sub("V", "", row_id)), year = observed_data$year[row_id])

ggplot(df_ppv_long, aes(factor(year), ppv)) +
geom_boxplot(outlier.shape = NA) +
scale_y_continuous(limits = c(0, 1), expand = c(0, 0), labels = scales::percent_format(accuracy = 1)) +
labs(x = "Year", y = "Positive predictive value") +
theme(axis.title = element_text(size = 17, family = "sans", colour = "black"),
      axis.text = element_text(size = 15, family = "sans", color = "black"))