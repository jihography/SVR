###############################################################################
## result_placebo.R
## ---------------------------------------------------------------------------
## IN-SPACE PLACEBO (Abadie permutation test)
##   각 control 을 가짜-처치로 두고 나머지 control 로 적합 →
##   post/pre RMSPE ratio 분포와 실제 treated band 비교.
##   treated/control 에 "동일한 단일-타깃 절차"를 적용해야 순열검정이 성립.
##
##   method 선택 (placebo_method):
##     "SMAC" : headline 과 동일 모형을 bands=1 로 재적합. 단, 단위가 1개라
##              treated_radius 공간 borrowing 은 작동하지 않아 사실상 단일
##              베이지안 vertical regression 으로 붕괴(= BVR 과 유사). 가장 비쌈.
##     "BVR"  : 단일 단위 베이지안 vertical regression (자연스러운 단일모드).
##     "SR"   : ridge. control 이 적고(여기선 8개) 선형종속이라 OLS 는 완전보간
##              (pre_RMSPE≈0)되어 ratio 퇴화 → 정규화로 방어. 빠르고 안정.
##     "OLS"  : sepOLS 와 동일(donor 충분할 때만).
##
##   입력: result_fit.R 가 저장한 res (+ t0, bands)
##   주의: control 8개 → permutation p 의 최소값은 1/(8+1)=0.111.
###############################################################################

setwd("/Users/kjh/Documents/innocity/SVR")

library(tidyverse)
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")

# ----------------------------- 설정/로드 ------------------------------------ #
# VARIANT 후보 (데이터 파일명 접미사 = unit_selection 의 SCENARIO + SELECTION_TAG):
#   "no_gumi_manual8greedy"  : 구미 제외 + 수동 8개 center (greedy 합산)  ← 현재
#   "no_gumi_densmatch"      : 구미 제외 + 밀도 percentile 자동매칭 control
#   "no_gumi"                : 구미 제외 (구버전)
#   "with_gumi"              : 구미 포함
#   "osrm"                   : OSRM 도로거리 기반 밴드
VARIANT <- "no_gumi_manual8greedy"                # result_fit.R 와 동일하게
SUFFIX  <- if (VARIANT == "with_gumi") "" else paste0("_", VARIANT)
RES_PATH <- sprintf("Output/apr_sims/Results/res_NA_t2013%s.RData", SUFFIX)
# 출력 디렉토리 (이미지/CSV) — 없으면 생성
OUT_DIR <- "Output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
outf <- function(name) file.path(OUT_DIR, sub("(\\.[A-Za-z0-9]+)$", paste0(SUFFIX, "\\1"), name))

load(RES_PATH)                                      # res, (t0, bands)
stopifnot(exists("res"))
if (!exists("t0") || !exists("bands")) {            # 구버전 RData 보강
  .env <- new.env()
  load(sprintf("SVR_input_data_t2013_%s.RData", VARIANT), envir = .env)
  if (!exists("t0"))    t0    <- .env$t0
  if (!exists("bands")) bands <- if (!is.null(.env$bands)) .env$bands else 5L
}
stopifnot(exists("t0"), exists("bands"))

years         <- as.integer(sub("^Y", "", rownames(res$sim)))
treated_names <- grep("Treated", colnames(res$sim), value = TRUE)
control_names <- setdiff(colnames(res$sim), treated_names)
n_control     <- length(control_names)

# ----------------------------- placebo 설정 --------------------------------- #
placebo_method <- "SMAC"    # "SMAC" / "BVR" / "SR"(ridge) / "OLS"
placebo_iter   <- 2000      # placebo Stan 표본 (point 예측용이라 본 적합보다 작게)
placebo_warm   <- 1000
placebo_chains <- 2
placebo_h      <- 0         # 단일 단위 spatial radius (treated/control 동일 적용)

# 베이지안 placebo(SMAC/BVR)는 rstan + 메서드 source + 컴파일 캐시 필요
if (placebo_method %in% c("SMAC", "BVR")) {
  suppressMessages(library(rstan))
  rstan_options(auto_write = TRUE)              # 컴파일 1회 후 재사용 (다회 fit 가속)
  iter <- placebo_iter; warm <- placebo_warm    # SMAC()/sepBVR() 가 global 로 읽음
  if (placebo_method == "SMAC") source("Methods/SMAC_method.R")
  if (placebo_method == "BVR")  source("Methods/sepBVR_method.R")
}

# 단일 타깃을 donor 들로 적합 → 표준화 잔차(gap) 기반 pre/post RMSPE
.placebo_fit_unit <- function(target, donors, t0, method = "OLS") {
  Tt  <- length(target)
  pre <- seq_len(t0); post <- (t0 + 1):Tt

  z  <- function(v, ref) (v - mean(ref)) / sd(ref)         # pre-기간 평균·SD 표준화
  tg <- z(target, target[pre])
  Dz <- apply(donors, 2, function(col) z(col, col[pre]))
  Dz[!is.finite(Dz)] <- 0                                  # 분산 0 donor 방어

  pred <- switch(method,
    "OLS" = {
      beta <- MASS::ginv(Dz[pre, , drop = FALSE]) %*% tg[pre]  # 절편 없음 (sepOLS 동일)
      as.numeric(Dz %*% beta)
    },
    "SR" = {
      suppressWarnings({                                   # 소표본 cv.glmnet 경고 억제
        cvg <- glmnet::cv.glmnet(Dz[pre, , drop = FALSE], tg[pre],
                                 alpha = 0, nfolds = min(10, length(pre)))
        bm  <- glmnet::glmnet(Dz[pre, , drop = FALSE], tg[pre], alpha = 0, lambda = cvg$lambda.min)
      })
      as.numeric(cbind(1, Dz) %*% as.matrix(coef(bm)))
    },
    "SMAC" = {                                             # 단일밴드 SMAC (headline 동일 모형)
      num_controls <<- ncol(Dz)                            # SMAC() 가 global 로 읽음
      tryCatch({
        fit <- SMAC(ym.pre = matrix(tg[pre], ncol = 1),
                    x.pre = Dz[pre, , drop = FALSE], x = Dz,
                    treated_radius = as.array(placebo_h),  # n_t=1 → 길이-1 벡터 강제
                    chains = placebo_chains)
        apply(fit$ynn[, , 1], 2, median)                   # ynn: draws x T x 1 → time별 median
      }, error = function(e) {
        message("  [SMAC] fit 실패 → NA: ", conditionMessage(e)); rep(NA_real_, Tt)
      })
    },
    "BVR" = {                                              # 단일 단위 베이지안 vertical reg
      tryCatch({
        bvr <- sepBVR(ym.pre = matrix(tg[pre], ncol = 1),
                      x.pre = Dz[pre, , drop = FALSE], x = Dz, chains = placebo_chains)
        apply(bvr[["treated 1"]]$y_new, 2, median)         # y_new: draws x T → time별 median
      }, error = function(e) {
        message("  [BVR] fit 실패 → NA: ", conditionMessage(e)); rep(NA_real_, Tt)
      })
    },
    stop("unknown placebo_method: ", method)
  )
  gap <- tg - pred                                          # 표준화 스케일 (ratio 는 스케일 불변)
  list(pre_RMSPE = sqrt(mean(gap[pre]^2)),
       post_RMSPE = sqrt(mean(gap[post]^2)), gap = gap)
}

set.seed(1)                                          # ridge cv.glmnet 폴드 재현성
sim_mat  <- as.matrix(res$sim)
ctrl_mat <- sim_mat[, control_names, drop = FALSE]

# treated: donor = 전체 control / control: donor = 나머지 control (leave-self-out)
fit_one <- function(nm, type) {
  message(sprintf("[placebo:%s] %s (%s)", placebo_method, nm, type))
  donors <- if (type == "treated") ctrl_mat
            else ctrl_mat[, setdiff(control_names, nm), drop = FALSE]
  f <- .placebo_fit_unit(sim_mat[, nm], donors, t0, placebo_method)
  f$unit <- nm; f$type <- type; f
}
fits <- c(lapply(treated_names, fit_one, type = "treated"),
          lapply(control_names, fit_one, type = "control"))

placebo_tbl <- map_dfr(fits, function(f) tibble(
  unit = f$unit, type = f$type,
  pre_RMSPE = f$pre_RMSPE, post_RMSPE = f$post_RMSPE,
  ratio = f$post_RMSPE / f$pre_RMSPE))

gap_long <- map_dfr(fits, function(f) tibble(
  unit = f$unit, type = f$type, year = years, gap = f$gap))

# permutation p-value: treated ratio 가 placebo(=control) 분포에서 차지하는 위치
ctrl_ratio <- placebo_tbl$ratio[placebo_tbl$type == "control"]
ctrl_ratio <- ctrl_ratio[is.finite(ctrl_ratio)]     # 실패한 fit(NA) 제외
treated_pl <- placebo_tbl %>%
  filter(type == "treated") %>%
  mutate(p_value = vapply(ratio, function(r) (1 + sum(ctrl_ratio >= r)) / (1 + length(ctrl_ratio)), numeric(1)))

placebo_out <- bind_rows(treated_pl,
                         placebo_tbl %>% filter(type == "control") %>% mutate(p_value = NA_real_)) %>%
  arrange(type, desc(ratio))

cat("\n== In-space placebo: treated band 별 RMSPE ratio & permutation p ==\n")
print(treated_pl)
write.csv(placebo_out, outf("placebo_inspace.csv"), row.names = FALSE)

# --- (A) gap 궤적: treated(강조) vs control(회색) ---
split_x <- years[t0] + 0.5
p_gap <- ggplot() +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = split_x, linetype = "dashed") +
  geom_line(data = filter(gap_long, type == "control"),
            aes(year, gap, group = unit), color = "gray75", linewidth = 0.4, alpha = 0.7) +
  geom_line(data = filter(gap_long, type == "treated"),
            aes(year, gap, color = unit), linewidth = 1.0) +
  scale_color_brewer(palette = "Set1", name = "Treated band",
                     labels = function(x) sub("Treated_", "", x)) +
  labs(x = "Year", y = "Gap (standardized residual)",
       title = "In-space placebo: treated vs control gaps",
       subtitle = sprintf("donor pool = %d controls | method = %s", n_control, placebo_method)) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

ggsave(outf("placebo_gaps.png"), p_gap, width = 11, height = 5, dpi = 500)

# --- (B) RMSPE ratio 순열 분포 + treated 표시 ---
p_ratio <- placebo_tbl %>%
  mutate(type = factor(type, levels = c("control", "treated"))) %>%
  ggplot(aes(x = ratio, y = type, color = type)) +
  geom_jitter(height = 0.12, width = 0, size = 2.6, alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = treated_pl,
    aes(x = ratio, y = "treated",
        label = sprintf("%s (p=%.3f)", sub("Treated_", "", unit), p_value)),
    inherit.aes = FALSE, size = 3, direction = "y", nudge_y = 0.3, segment.alpha = 0.4) +
  scale_color_manual(values = c(control = "gray60", treated = "#009E73"), guide = "none") +
  labs(x = "post/pre RMSPE ratio", y = NULL,
       title = "Permutation distribution of RMSPE ratio",
       subtitle = "treated band 가 control placebo 분포의 오른쪽 꼬리에 있을수록 효과 유의") +
  theme_bw(base_size = 14)

ggsave(outf("placebo_ratio.png"), p_ratio, width = 11, height = 4, dpi = 500)

message("placebo 완료: placebo_inspace.csv, placebo_gaps.png, placebo_ratio.png (method=",
        placebo_method, ", SUFFIX='", SUFFIX, "')")
