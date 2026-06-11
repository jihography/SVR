###############################################################################
## result_postestimation.R
## ---------------------------------------------------------------------------
## result_fit.R 가 저장한 `res` 를 불러, 후처리 일체를 수행한다:
##   1) beta 비교 그림 (Figure 4) + 변동(SD) 요약 + boxplot
##   2) treatment-effect 그림 (Figure 5) + intertemporal average
##   3) 모형 비교 통계 (fit_summary: 행=지표, 열=모델)
##   4) poster 시각화 2종: (A) Core 관측 vs 합성통제  (B) 밴드 간 가중치 상관
## placebo 는 result_placebo.R, DiD 는 result_did.R 로 분리.
###############################################################################

setwd("/Users/kjh/Documents/innocity/SVR")

library(tidyverse)
library(patchwork)
library(purrr)
library(tibble)
library(tidyr)
if (!requireNamespace("ggh4x", quietly = TRUE)) install.packages("ggh4x")
library(ggh4x)

# ----------------------------- 설정/로드 ------------------------------------ #
# VARIANT 후보 (데이터 파일명 접미사 = unit_selection 의 SCENARIO + SELECTION_TAG):
#   "no_gumi_manual8greedy"  : 구미 제외 + 수동 8개 center (greedy 합산)  ← 현재
#   "no_gumi_densmatch"      : 구미 제외 + 밀도 percentile 자동매칭 control
#   "no_gumi"                : 구미 제외 (구버전)
#   "with_gumi"              : 구미 포함
#   "osrm"                   : OSRM 도로거리 기반 밴드
VARIANT <- "no_gumi"                  # result_fit.R 와 동일하게
SUFFIX  <- if (VARIANT == "with_gumi") "" else paste0("_", VARIANT)
RES_PATH <- sprintf("Output/apr_sims/Results/res_NA_t2013%s.RData", SUFFIX)

# 출력 디렉토리 (이미지/CSV) — 없으면 생성
OUT_DIR <- "Output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 출력 파일명에 VARIANT 접미사 삽입 (확장자 앞) + OUT_DIR 하위에 저장
outf <- function(name) file.path(OUT_DIR, sub("(\\.[A-Za-z0-9]+)$", paste0(SUFFIX, "\\1"), name))

load(RES_PATH)                                      # res, (t0, bands)
stopifnot(exists("res"))

# 구 버전 RData 는 res 만 저장했을 수 있음 → 입력 데이터에서 t0/bands 보강
if (!exists("t0") || !exists("bands")) {
  .env <- new.env()
  load(sprintf("SVR_input_data_t2013_%s.RData", VARIANT), envir = .env)
  if (!exists("t0"))    t0    <- .env$t0
  if (!exists("bands")) bands <- if (!is.null(.env$bands)) .env$bands else 5L
}
stopifnot(exists("t0"), exists("bands"))

# 공통 파생값
years         <- as.integer(sub("^Y", "", rownames(res$sim)))
treated_names <- grep("Treated", colnames(res$sim), value = TRUE)
control_names <- setdiff(colnames(res$sim), treated_names)
n_treated     <- length(treated_names)
n_control     <- length(control_names)
stopifnot(n_treated > 0, n_control > 0)

COL <- c(BVR = "#F0E442", OLS = "#56B4E9", SVR = "#009E73")

band_levels <- c("Treated_Core", "Treated_Band_5km", "Treated_Band_10km",
                 "Treated_Band_15km", "Treated_Band_20km")


## ===========================================================================
## 1) BETA 비교 (Figure 4)
## ===========================================================================

# BVR: treated 1..K 의 gp1(draw x Kcol)에서 control 별 posterior median 행렬 만들기
get_bvr_beta_matrix_from_gp1 <- function(res, n_treated, treated_names, control_names = NULL) {
  BVR <- res$est$BVR
  treated_slots <- grep("^treated\\s*[0-9]+$", names(BVR), value = TRUE)
  if (length(treated_slots) == 0) stop("BVR에서 'treated 1' 형태 키를 못 찾음.")
  ord <- order(as.integer(gsub("^treated\\s*", "", treated_slots)))
  treated_slots <- treated_slots[ord]

  gp1_1 <- as.matrix(BVR[[treated_slots[1]]]$gp1)   # draws x K
  K <- ncol(gp1_1)

  drop_idx <- integer(0)
  if (!is.null(colnames(gp1_1))) {
    drop_idx <- which(tolower(colnames(gp1_1)) %in% c("(intercept)", "intercept", "beta0", "const"))
  }
  if (length(drop_idx) == 0 && !is.null(control_names) && K == length(control_names) + 1) {
    drop_idx <- 1
  }

  gp1_1_use <- gp1_1
  if (length(drop_idx) > 0) gp1_1_use <- gp1_1[, -drop_idx, drop = FALSE]
  K2 <- ncol(gp1_1_use)

  T_use <- min(n_treated, length(treated_slots))
  B <- matrix(NA_real_, nrow = T_use, ncol = K2)
  for (t in seq_len(T_use)) {
    gp1 <- as.matrix(BVR[[treated_slots[t]]]$gp1)
    if (length(drop_idx) > 0) gp1 <- gp1[, -drop_idx, drop = FALSE]
    B[t, ] <- apply(gp1, 2, median)
  }
  rownames(B) <- treated_names[seq_len(nrow(B))]
  if (!is.null(colnames(gp1_1_use))) {
    colnames(B) <- colnames(gp1_1_use)
  } else if (!is.null(control_names) && length(control_names) == ncol(B)) {
    colnames(B) <- control_names
  } else {
    colnames(B) <- paste0("C", seq_len(ncol(B)))
  }
  B
}

B_bvr <- get_bvr_beta_matrix_from_gp1(res, n_treated, treated_names, control_names)

# OLS: (control x treated) -> (treated x control)
B_ols <- t(as.matrix(res$est$OLS))
rownames(B_ols) <- treated_names[seq_len(nrow(B_ols))]
colnames(B_ols) <- control_names[seq_len(ncol(B_ols))]

# SVR(SMAC): draw x treated x control -> posterior median (treated x control)
beta_draws <- res$est$SMAC$beta
stopifnot(length(dim(beta_draws)) == 3, dim(beta_draws)[2] == n_treated)
B_svr <- apply(beta_draws, c(2, 3), median)
rownames(B_svr) <- treated_names
if (is.null(colnames(B_svr)) && length(control_names) == ncol(B_svr)) {
  colnames(B_svr) <- control_names
}

# 세 방법의 control 차원 일치
if (!is.null(colnames(B_bvr)) && !is.null(colnames(B_ols)) && !is.null(colnames(B_svr))) {
  common_controls <- Reduce(intersect, list(colnames(B_bvr), colnames(B_ols), colnames(B_svr)))
  if (length(common_controls) == 0) stop("세 방법의 control colnames 교집합이 0.")
  B_bvr <- B_bvr[, common_controls, drop = FALSE]
  B_ols <- B_ols[, common_controls, drop = FALSE]
  B_svr <- B_svr[, common_controls, drop = FALSE]
} else {
  k <- min(ncol(B_bvr), ncol(B_ols), ncol(B_svr))
  B_bvr <- B_bvr[, 1:k, drop = FALSE]
  B_ols <- B_ols[, 1:k, drop = FALSE]
  B_svr <- B_svr[, 1:k, drop = FALSE]
}

beta_long <- bind_rows(
  as_tibble(as.data.frame(as.table(B_bvr))) %>% mutate(method = "BVR"),
  as_tibble(as.data.frame(as.table(B_ols))) %>% mutate(method = "OLS"),
  as_tibble(as.data.frame(as.table(B_svr))) %>% mutate(method = "SVR")
) %>%
  rename(treated = Var1, control = Var2, beta = Freq) %>%
  mutate(method = factor(method, levels = c("BVR", "OLS", "SVR")))

# 패널 제목(facet strip)은 control(대조지역) 변수명. 한글로 바꾸려면
# 아래 named vector에 "원본이름" = "한글이름" 매핑을 채워 넣으면 된다.
control_labs_ko <- c(
  # "control_name_1" = "대조지역1",
)
control_labeller <- if (length(control_labs_ko) > 0) {
  labeller(control = control_labs_ko)
} else "label_value"

p_fig4 <- beta_long %>%
  ggplot(aes(x = treated, y = beta, color = method, group = method)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ control, nrow = 2, labeller = control_labeller) +
  scale_color_manual(values = COL, name = "방법") +
  scale_x_discrete(breaks = treated_names, labels = as.character(seq_along(treated_names))) +
  labs(x = "처치 지역", y = expression(beta[c])) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(outf("models_2013.png"), p_fig4, width = 13, height = 6, dpi = 500)

# --- beta 변동 point plot ---
beta_summary <- beta_long %>%
  group_by(method, control) %>%
  summarise(mean = mean(beta, na.rm = TRUE), sd = sd(beta, na.rm = TRUE),
            min = min(beta, na.rm = TRUE),   max = max(beta, na.rm = TRUE),
            .groups = "drop")

dodge <- position_dodge(width = 0.65)
p_fig4_pts <- ggplot() +
  geom_errorbar(data = beta_summary,
                aes(x = control, ymin = min, ymax = max, color = method),
                position = dodge, width = 0.35, linewidth = 0.7, alpha = 0.5) +
  geom_point(data = beta_long, aes(x = control, y = beta, color = method),
             position = dodge, size = 2.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = COL, name = "Method") +
  labs(x = NULL, y = expression(beta[c]),
       title = "Per-control beta across 5 treated bands",
       subtitle = "Range bar = min-max | tight cluster = stable estimates") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(outf("models_pts.png"), p_fig4_pts, width = 14, height = 6, dpi = 500)

# --- beta 변동(SD) 숫자 요약 ---
beta_sd_summary <- beta_summary %>%
  group_by(method) %>%
  summarise(mean_sd_across_controls = mean(sd, na.rm = TRUE),
            max_sd = max(sd, na.rm = TRUE),
            median_range = median(max - min, na.rm = TRUE)) %>%
  arrange(mean_sd_across_controls)

cat("== beta 변동 종합 (작을수록 안정) ==\n"); print(beta_sd_summary)
write.csv(beta_sd_summary, outf("beta_sd_summary.csv"), row.names = FALSE)

# --- pooled boxplot ---
p_box <- ggplot(beta_long, aes(x = method, y = beta, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = COL) +
  labs(x = NULL, y = expression(beta[c]),
       title = "Pooled distribution of beta (across all treated x control)") +
  theme_bw(base_size = 14) + theme(legend.position = "none")

ggsave(outf("models_box.png"), p_box, width = 6, height = 5, dpi = 500)


## ===========================================================================
## 2) TREATMENT EFFECT 그림 (Figure 5) + intertemporal average
## ===========================================================================
COL_SVR <- "#009E73"

y_obs <- as.matrix(res$sim)
if (ncol(y_obs) != n_treated) y_obs <- y_obs[, seq_len(n_treated), drop = FALSE]

eff_med <- -1 * as.matrix(res$point$bias$SMACmedian)   # point estimate
mu_lo   <- as.matrix(res$ci$SMAC$lower_bound)
mu_up   <- as.matrix(res$ci$SMAC$upper_bound)
eff_lo  <- y_obs - mu_lo
eff_up  <- y_obs - mu_up
stopifnot(all(dim(eff_med) == dim(y_obs)),
          all(dim(eff_lo)  == dim(y_obs)),
          all(dim(eff_up)  == dim(y_obs)))

te_df <- purrr::map_dfr(seq_along(band_levels), function(j) {
  nm <- band_levels[j]
  tibble(treated = nm, year = years,
         median = eff_med[, nm], lower = eff_lo[, nm], upper = eff_up[, nm])
}) %>%
  mutate(treated = factor(treated, levels = band_levels))

y_scales <- list(
  scale_y_continuous(limits = c(-5, 55)),    # Treated_Core
  scale_y_continuous(limits = c(-20, 20)),   # 5km
  scale_y_continuous(limits = c(-20, 20)),   # 10km
  scale_y_continuous(limits = c(-20, 20)),     # 15km
  scale_y_continuous(limits = c(-20, 20))      # 20km
)

# 패널 제목(facet strip) 한글 라벨
band_labs_ko <- c(
  "Treated_Core"      = "율곡동(경북혁신도시)",
  "Treated_Band_5km"  = "5km 밴드",
  "Treated_Band_10km" = "10km 밴드",
  "Treated_Band_15km" = "15km 밴드",
  "Treated_Band_20km" = "20km 밴드"
)

p_fig5 <- te_df %>%
  ggplot(aes(x = year, y = median)) +
  theme(aspect_ratio = 0.2) +
  scale_x_continuous(breaks = c(1997, 2007, 2013)) +
  geom_hline(yintercept = 0, color = "gray60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = COL_SVR, alpha = 0.25) +
  geom_line(color = COL_SVR, linewidth = 1.1) +
  geom_vline(xintercept = 2013, linewidth = 0.7) +
  facet_wrap(~ treated, nrow = 1, scales = "free_y",
             labeller = labeller(treated = band_labs_ko)) +
  ggh4x::facetted_pos_scales(y = y_scales) +
  labs(x = "연도", y = "처치 효과") +
  theme_bw(base_size = 14) + theme(panel.grid.minor = element_blank())

ggsave(outf("post_2013.png"), p_fig5, width = 14, height = 3.6, dpi = 500)

# --- intertemporal average (draw 별 사후 평균 효과 → median + CrI) ---
intertemporal_avg_from_beta <- function(res, t0, probs = c(0.025, 0.5, 0.975),
                                        treated_names = NULL, control_names = NULL) {
  sim <- as.matrix(res$sim); Tt <- nrow(sim)

  if (is.null(treated_names)) {
    cand <- band_levels
    treated_names <- if (all(cand %in% colnames(sim))) cand else grep("Treated", colnames(sim), value = TRUE)
  }
  I <- length(treated_names); stopifnot(I > 0)
  if (is.null(control_names)) control_names <- setdiff(colnames(sim), treated_names)

  beta_raw  <- as.array(res$est$SMAC$beta)
  beta0_raw <- as.array(res$est$SMAC$beta0)
  d <- dim(beta_raw); stopifnot(length(d) == 3)
  S <- max(d); I0 <- min(d); C0 <- setdiff(d, c(S, I0))[1]
  if (I0 != I) stop("treated 개수 불일치: beta=", I0, ", names=", I)

  idx_draw  <- which(d == S)[1]; idx_treat <- which(d == I)[1]
  idx_ctrl  <- setdiff(1:3, c(idx_draw, idx_treat))[1]
  beta <- aperm(beta_raw, c(idx_draw, idx_treat, idx_ctrl))   # [S,I,C]

  d0 <- dim(beta0_raw)
  if (length(d0) == 2) {
    idx_draw0 <- which(d0 == S)[1]; idx_treat0 <- which(d0 == I)[1]
    beta0 <- if (idx_draw0 == 1 && idx_treat0 == 2) beta0_raw else t(beta0_raw)
  } else if (length(d0) == 1 && length(beta0_raw) == S * I) {
    beta0 <- matrix(beta0_raw, nrow = S, ncol = I, byrow = TRUE)
  } else stop("beta0 dim 해석 불가: ", paste(d0, collapse = "x"))

  post_idx <- (t0 + 1):Tt; post_len <- length(post_idx); stopifnot(post_len > 0)
  Y_post <- sim[post_idx, treated_names, drop = FALSE]

  if (length(control_names) < C0) stop("control_names(", length(control_names), ") < beta controls(", C0, ")")
  if (length(control_names) > C0) control_names <- control_names[seq_len(C0)]
  X_post <- sim[post_idx, control_names, drop = FALSE]

  Delta_draws <- matrix(NA_real_, nrow = S, ncol = I)
  for (i in seq_len(I)) {
    B_i  <- beta[, i, , drop = FALSE][, 1, ]   # [S,C]
    pred <- X_post %*% t(B_i)                  # [post,S]
    pred <- sweep(pred, 2, beta0[, i], "+")
    tau  <- matrix(Y_post[, i], nrow = post_len, ncol = S) - pred
    Delta_draws[, i] <- colMeans(tau)
  }

  tibble(band = treated_names, post_len = post_len,
         mean = colMeans(Delta_draws), median = apply(Delta_draws, 2, median),
         lower = apply(Delta_draws, 2, quantile, probs[1]),
         upper = apply(Delta_draws, 2, quantile, probs[3]))
}

avg_df <- intertemporal_avg_from_beta(res, t0 = t0)
print(avg_df)
write.csv(avg_df,   outf("post_t2013.csv"), row.names = FALSE)
write.csv(res$sim,  outf("sim.csv"),        row.names = TRUE)   # DiD 입력


## ===========================================================================
## 3) 모형 비교 통계 (fit_summary: 행=지표, 열=모델)
##    post_estimate / post_CI_width 는 Core + 5km band 만 사용.
## ===========================================================================
Tt       <- nrow(res$sim)
Y_full   <- res$sim[, 1:bands, drop = FALSE]
Y_pre    <- Y_full[1:t0, , drop = FALSE]
Y_post   <- Y_full[(t0 + 1):Tt, , drop = FALSE]
post_idx <- (t0 + 1):Tt

preds <- list(SC = res$cal$SC, SR = res$cal$SR, OLS = res$cal$OLS,
              BVR = res$cal$BVRmedian, BSC = res$cal$BSCmedian, SMAC = res$cal$SMACmedian)

rmse_per_band <- function(Y, Yhat) sqrt(colMeans((Y - Yhat)^2, na.rm = TRUE))
mae_per_band  <- function(Y, Yhat) colMeans(abs(Y - Yhat), na.rm = TRUE)
r2_per_band   <- function(Y, Yhat) {
  ss_res <- colSums((Y - Yhat)^2, na.rm = TRUE)
  ss_tot <- colSums(sweep(Y, 2, colMeans(Y, na.rm = TRUE))^2, na.rm = TRUE)
  1 - ss_res / ss_tot
}

# (a) band 단위 fit 통계
fit_stats <- map_dfr(names(preds), function(m) {
  Yhat    <- preds[[m]]
  te_post <- Y_post - Yhat[post_idx, , drop = FALSE]
  tibble(
    method        = m,
    band          = colnames(Y_full),
    pre_RMSPE     = rmse_per_band(Y_pre,  Yhat[1:t0, , drop = FALSE]),
    post_RMSPE    = rmse_per_band(Y_post, Yhat[post_idx, , drop = FALSE]),
    ratio         = post_RMSPE / pre_RMSPE,
    pre_MAE       = mae_per_band( Y_pre,  Yhat[1:t0, , drop = FALSE]),
    pre_R2        = r2_per_band(  Y_pre,  Yhat[1:t0, , drop = FALSE]),
    post_estimate = colMeans(te_post, na.rm = TRUE)
  )
})

# (b) 추론 정밀도: CI 폭 + pre-coverage. post_CI_width 는 Core + 5km 만.
core_bands <- c("Treated_Core", "Treated_Band_5km")

inference_stats <- map_dfr(names(res$ci), function(m) {
  ci_m  <- res$ci[[m]];  if (is.null(ci_m)) return(NULL)
  cov_m <- res$coverage[[m]]
  lb <- as.matrix(ci_m$lower_bound); ub <- as.matrix(ci_m$upper_bound)
  width <- ub - lb

  core_cols  <- intersect(core_bands, colnames(width))
  width_core <- if (length(core_cols) > 0) width[, core_cols, drop = FALSE] else width

  pre_w  <- mean(width[1:t0, ],          na.rm = TRUE)
  post_w <- mean(width_core[post_idx, ], na.rm = TRUE)

  if (isTRUE(all(width == 0, na.rm = TRUE))) {     # OLS 처럼 점추정만 → NA
    pre_w <- NA_real_; post_w <- NA_real_; valid_ci <- FALSE
  } else valid_ci <- TRUE

  tibble(method = m, pre_ci_width = pre_w, post_CI_width = post_w,
         pre_coverage = if (valid_ci && !is.null(cov_m))
                          mean(as.matrix(cov_m)[1:t0, ], na.rm = TRUE) else NA_real_)
})

# 통합 (모델별 1행) → post_estimate 는 Core + 5km 평균
fit_summary_long <- fit_stats %>%
  group_by(method) %>%
  summarise(
    pre_RMSPE     = mean(pre_RMSPE,  na.rm = TRUE),
    pre_R2        = mean(pre_R2,     na.rm = TRUE),
    post_estimate = mean(post_estimate[band %in% core_bands], na.rm = TRUE),
    post_RMSPE    = mean(post_RMSPE, na.rm = TRUE),
    ratio         = mean(ratio,      na.rm = TRUE),
    pre_MAE       = mean(pre_MAE,    na.rm = TRUE)
  ) %>%
  left_join(inference_stats, by = "method") %>%
  select(method, pre_RMSPE, pre_R2, post_estimate, post_CI_width,
         post_RMSPE, ratio, pre_MAE, pre_ci_width, pre_coverage)

# 전치: 행=지표, 열=모델
fit_summary <- fit_summary_long %>%
  pivot_longer(-method, names_to = "stat") %>%
  pivot_wider(names_from = method, values_from = value) %>%
  slice(match(c("pre_RMSPE", "pre_R2", "post_estimate", "post_CI_width",
                "post_RMSPE", "ratio", "pre_MAE", "pre_ci_width", "pre_coverage"), stat))

# 보조 표: band x method pre-RMSPE
fit_wide_rmspe <- fit_stats %>%
  select(method, band, pre_RMSPE) %>%
  pivot_wider(names_from = method, values_from = pre_RMSPE)

cat("\n== Pre-treatment RMSPE (band x method) ==\n"); print(fit_wide_rmspe)
cat("\n== Method summary (행=지표, 열=모델) ==\n");    print(fit_summary)

write.csv(fit_stats,   outf("fit_stats.csv"),   row.names = FALSE)
write.csv(fit_summary, outf("fit_summary.csv"), row.names = FALSE)



## ===========================================================================
## 4) POSTER 시각화 2종
##    (A) 율곡동(Core): 관측 vs 합성통제 fit
##    (B) SVR(SMAC)이 추정한 밴드 간 가중치 공간 상관 (5×5 히트맵)
## ===========================================================================
COL_NEG    <- "#C0584B"   # (B) 음의 상관용
band_short <- c("Core", "5km", "10km", "15km", "20km")

## ------------------------------------------------------------
## (A) Synthetic fit checking (Abadie 스타일): 관측 vs 합성통제 — 원레벨(인구/값)
##     SVR 적합(res$cal)은 detrended 스케일이라 DiD gap 처럼 보임 → 실제 레벨로
##     역변환해 그린다. detrend 정의 (unit_selection.Rmd):
##        sim = (raw - pre_mean_i)/pre_sd_i - trend_t ,  trend_t = mean(표준화 control)
##     역변환:  raw           = (sim + trend_t)*pre_sd_i + pre_mean_i
##              synthetic_raw = (cal + trend_t)*pre_sd_i + pre_mean_i   (cal = res$cal$SMACmedian)
## ------------------------------------------------------------
TREAT_YEAR <- 2013

# 원레벨(sim_orig) + detrend 파라미터 복원
.ie <- new.env()
load(sprintf("SVR_input_data_t2013_%s.RData", VARIANT), envir = .ie)
sim_orig  <- as.matrix(.ie$sim_orig)
pre_idx   <- which(years < TREAT_YEAR)
pre_mean  <- colMeans(sim_orig[pre_idx, , drop = FALSE], na.rm = TRUE)
pre_sd    <- apply(sim_orig[pre_idx, , drop = FALSE], 2, sd, na.rm = TRUE)
Y_std     <- sweep(sweep(sim_orig, 2, pre_mean, "-"), 2, pre_sd, "/")
ctrl_cols <- setdiff(colnames(sim_orig), treated_names)
trend_t   <- rowMeans(Y_std[, ctrl_cols, drop = FALSE], na.rm = TRUE)

cal_mat <- as.matrix(res$cal$SMACmedian)
colnames(cal_mat) <- colnames(res$sim)[seq_len(bands)]      # cal 은 colname 없음 → 부여

# 밴드별 관측/합성통제 (원레벨로 역변환)
fit_long <- purrr::map_dfr(band_levels, function(nm) {
  tibble(band = nm, year = years,
         Observed  = sim_orig[, nm],
         Synthetic = (cal_mat[, nm] + trend_t) * pre_sd[nm] + pre_mean[nm])
}) %>%
  pivot_longer(c(Observed, Synthetic), names_to = "series", values_to = "value") %>%
  mutate(band   = factor(band, levels = band_levels),
         series = factor(series, levels = c("Observed", "Synthetic")))

xbreaks <- c(1997, 2007, 2013)   # post_2013 플롯과 동일하게 일부 연도만 표시

# 관측=검정 실선 / 합성통제=초록 긴점선 (Abadie 구도 + color)
abadie_fit <- function(df) {
  ggplot(df, aes(year, value, color = series, linetype = series)) +
    geom_vline(xintercept = TREAT_YEAR, linetype = "dotted", linewidth = 0.6) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = c(Observed = "#1A1A1A", Synthetic = "#009E73"), name = NULL) +
    scale_linetype_manual(values = c(Observed = "solid", Synthetic = "longdash"), name = NULL) +
    scale_x_continuous(breaks = xbreaks) +
    scale_y_continuous(labels = scales::comma,
                       breaks = scales::pretty_breaks(n = 5)) +  # y축 간격 균등하게
    labs(x = "Year", y = "Outcome (level)") +
    theme_classic(base_size = 14)
}

# (A-1) Core(율곡동) 단일 패널 — 첨부 이미지와 동일 구도 (처치선 + 화살표 주석)
core_df <- filter(fit_long, band == "Treated_Core")
.ylo <- min(core_df$value, na.rm = TRUE) + 0.12 * diff(range(core_df$value, na.rm = TRUE))
p_core <- abadie_fit(core_df) +
  annotate("segment", x = TREAT_YEAR - 6, xend = TREAT_YEAR - 0.4, y = .ylo, yend = .ylo,
           arrow = grid::arrow(length = grid::unit(0.18, "cm")), linewidth = 0.5, color = "grey30") +
  annotate("text", x = TREAT_YEAR - 6, y = .ylo, label = "처치 2013",
           hjust = 1.1, size = 3.7, color = "grey30") +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(color = "black", linewidth = 0.3))

ggsave(outf("synthetic_fit_core.png"), p_core, width = 6.5, height = 4.5, dpi = 500)

# (A-2) 5개 밴드 전체 fit checking (facet, pre-기간 적합도 점검)
p_bands <- abadie_fit(fit_long) +
  facet_wrap(~ band, nrow = 1, scales = "free_y",
             labeller = as_labeller(setNames(band_short, band_levels))) +
  theme(legend.position = "top")

ggsave(outf("synthetic_fit_bands.png"), p_bands, width = 14, height = 3.6, dpi = 500)

## ------------------------------------------------------------
## (B) 밴드 간 가중치 공간 상관 — SVR이 반영하는 공간 의존성
##     B_svr: treated(5) × control 가중치 (posterior median)
##     cor(t(B_svr)) = 밴드끼리 통제군 가중치 벡터가 얼마나 닮았는가
## ------------------------------------------------------------
B_svr <- apply(res$est$SMAC$beta, c(2, 3), median)          # [treated, control]
rownames(B_svr) <- treated_names[seq_len(nrow(B_svr))]
B_svr <- B_svr[band_levels, , drop = FALSE]                 # 밴드 순서 정렬

Cmat <- cor(t(B_svr))                                       # 5×5
dimnames(Cmat) <- list(band_short, band_short)

cdf <- as.data.frame(as.table(Cmat)) %>%
  rename(b1 = Var1, b2 = Var2, corr = Freq) %>%
  mutate(b1 = factor(b1, levels = band_short),
         b2 = factor(b2, levels = rev(band_short)))         # y축 위→아래 Core..20km

p_dep <- ggplot(cdf, aes(b1, b2, fill = corr)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.2f", corr)), size = 4.3, color = "grey15") +
  scale_fill_gradient2(low = COL_NEG, mid = "white", high = COL_SVR,
                       midpoint = 0, limits = c(-1, 1), name = "상관") +
  coord_equal() +
  labs(x = NULL, y = NULL,
       title = "밴드 간 가중치 공간 상관",
       subtitle = "인접 밴드일수록 강한 의존 → SVR이 데이터로부터 추정") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        plot.subtitle = element_text(color = "grey35", size = 11))

ggsave(outf("spatial_dependence_corr.png"), p_dep, width = 5.4, height = 5.6, dpi = 500)

message("postestimation 완료: 그림/표 저장됨 (SUFFIX='", SUFFIX, "')")

