###############################################################################
## result_fit.R
## ---------------------------------------------------------------------------
## 6개 방법(SC, SR, OLS, BVR, BSC, SMAC) 적합 후, 중요 결과인 `res` 만 RData 로
## 저장한다. (구) result_orig.Rmd 의 "analysis" 청크를 분리한 것.
## 무거운 적합은 여기서만 1회 수행하고, 이후 plot/통계/placebo/DiD 는
## result_postestimation.R, result_did.R 에서 이 RData 를 불러 사용한다.
###############################################################################

setwd("/Users/kjh/Documents/innocity/SVR")

# ----------------------------- 설정 ----------------------------------------- #
# VARIANT 후보 (= unit_selection 의 SCENARIO + SELECTION_TAG, 데이터 파일명 접미사):
#   "no_gumi_manual8greedy"  : 구미 제외 + 수동 8개 center (greedy 합산)  ← 현재
#   "no_gumi_densmatch"      : 구미 제외 + 밀도 percentile 자동매칭 control
#   "no_gumi"                : 구미 제외 (구버전)
#   "with_gumi"              : 구미 포함
#   "osrm"                   : OSRM 도로거리 기반 밴드
VARIANT <- "no_gumi"                # unit_selection: SCENARIO + SELECTION_TAG
SUFFIX  <- if (VARIANT == "with_gumi") "" else paste0("_", VARIANT)
DATA_PATH <- sprintf("SVR_input_data_t2013_%s.RData", VARIANT)
RES_PATH  <- sprintf("Output/apr_sims/Results/res_NA_t2013%s.RData", SUFFIX)

# 추정 하이퍼파라미터
iter   <- 4000
warm   <- 2000
chains <- 4
treated_radius <- c(0, 0.25, 0.5, 0.75, 1)
method <- c("SC", "SR", "OLS", "BVR", "BSC", "SMAC")

# ------------------------- 함수 source -------------------------------------- #
source("Functions/sim_model_function.R")
source("Functions/estimation_function.R")
source("Functions/calculation_function.R")
source("Functions/point_estimate_function.R")
source("Functions/ci_function.R")
source("Functions/coverage_function.R")
source("Methods/helper/preset_function.R")
source("Methods/helper/ci_shen.R")
source("Methods/helper/ci_bayes.R")
source("Methods/helper/cover_helper.R")
source("Methods/helper/binder.R")
source("Methods/helper/wrapper.R")
source("Methods/helper/helper_vertical_regression.R")
source("Methods/helper/SCM_function.R")
source("Methods/sepSC_method.R")
source("Functions/helper/sepSC_calc_function.R")
source("Methods/sepSR_method.R")
source("Functions/helper/sepSR_calc_function.R")
source("Methods/sepOLS_method.R")
source("Methods/sepBVR_method.R")
source("Methods/sepBSC_method.R")
source("Methods/SMAC_method.R")

library(LowRankQP)
library(glmnet)
library(rstan)
library(spatstat)
library(Matrix)
library(fungible)
library(purrr)
rstan_options(auto_write = FALSE)

# --------------------------- 데이터 로드 ------------------------------------ #
external_env <- new.env()
loaded_vars  <- load(DATA_PATH, envir = external_env)
if (!"sim" %in% loaded_vars) stop("입력 파일에 'sim' 객체가 없습니다: ", DATA_PATH)

sim   <- as.matrix(external_env$sim)
bands <- if ("bands" %in% loaded_vars) external_env$bands else 5L
t0    <- if ("t0"    %in% loaded_vars) external_env$t0    else stop("t0 가 입력 파일에 없습니다.")
num_controls <- if ("num_controls" %in% loaded_vars) external_env$num_controls else ncol(sim) - bands
time_periods <- nrow(sim)            # ci_bayes() 가 전역으로 참조 (원본 Rmd 와 동일)

message(sprintf("Loaded sim: %d x %d | bands=%d | t0=%d | controls=%d",
                nrow(sim), ncol(sim), bands, t0, num_controls))

# 유효성 검사
stopifnot(bands >= 1, bands < ncol(sim), t0 >= 1, t0 < nrow(sim))
if (ncol(sim) - bands != num_controls) {
  num_controls <- ncol(sim) - bands
  message("num_controls 를 데이터에 맞춰 재설정: ", num_controls)
}

# ----------------------------- 적합 ---------------------------------------- #
est     <- estimation(sim = sim, t0 = t0, bands = bands, iter = iter, warm = warm,
                      norm = TRUE, method = method,
                      treated_radius = treated_radius, chains = chains)
cal     <- calculation(sim = sim, est = est, bands = bands, norm = TRUE)
point   <- point_estimate(sim, cal)                 # bias / squared error
c_interv<- ci(sim = sim, est = est, cal = cal, t0 = t0, norm = TRUE)
cover   <- coverage(sim, interv = c_interv)

# --------------------- 결과 저장 (res + 해석에 필요한 스칼라) ---------------- #
res <- list(sim = sim, est = est, cal = cal, beta_true = NULL,
            point = point, ci = c_interv, coverage = cover)

dir.create(dirname(RES_PATH), recursive = TRUE, showWarnings = FALSE)
save(res, t0, bands, file = RES_PATH)   # t0, bands 는 res 해석에 필수인 작은 스칼라
message("Saved fitted result -> ", RES_PATH)
