###############################################################################
## result_did.R
## ---------------------------------------------------------------------------
## TWFE event-study DiD (밴드별). 입력은 result_postestimation.R 가 저장한
## sim.csv (time x units). (구) result_orig.Rmd 의 "DiD" 청크를 분리한 것.
###############################################################################

setwd("/Users/kjh/Documents/innocity/SVR")

library(dplyr)
library(tidyr)
library(ggplot2)
library(fixest)
library(broom)

# ----------------------------- 설정 ----------------------------------------- #
# VARIANT 후보 (데이터 파일명 접미사 = unit_selection 의 SCENARIO + SELECTION_TAG):
#   "no_gumi_manual8greedy"  : 구미 제외 + 수동 8개 center (greedy 합산)  ← 현재
#   "no_gumi_densmatch"      : 구미 제외 + 밀도 percentile 자동매칭 control
#   "no_gumi"                : 구미 제외 (구버전)
#   "with_gumi"              : 구미 포함
#   "osrm"                   : OSRM 도로거리 기반 밴드
VARIANT <- "no_gumi_manual8greedy"                  # postestimation 과 동일하게
SUFFIX  <- if (VARIANT == "with_gumi") "" else paste0("_", VARIANT)
# 출력 디렉토리 (이미지/CSV) — postestimation 과 동일하게 Output/ 하위 사용
OUT_DIR <- "Output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
outf    <- function(name) file.path(OUT_DIR, sub("(\\.[A-Za-z0-9]+)$", paste0(SUFFIX, "\\1"), name))

TREAT_YEAR   <- 2013
treated_cols <- c("Treated_Core", "Treated_Band_5km", "Treated_Band_10km",
                  "Treated_Band_15km", "Treated_Band_20km")

# --------------------------- sim 불러오기 ----------------------------------- #
sim_df <- read.csv(outf("sim.csv"), check.names = FALSE)
rownames(sim_df) <- sim_df[[1]]
sim_df[[1]] <- NULL
years <- as.integer(sub("^Y", "", rownames(sim_df)))
stopifnot(all(treated_cols %in% colnames(sim_df)))

# --------------------------- long panel ------------------------------------- #
did_df <- sim_df %>%
  mutate(year = years) %>%
  pivot_longer(-year, names_to = "unit_id", values_to = "outcome") %>%
  mutate(
    rel_time = year - TREAT_YEAR,
    core = as.integer(unit_id == "Treated_Core"),
    b5   = as.integer(unit_id == "Treated_Band_5km"),
    b10  = as.integer(unit_id == "Treated_Band_10km"),
    b15  = as.integer(unit_id == "Treated_Band_15km"),
    b20  = as.integer(unit_id == "Treated_Band_20km"),
    band = case_when(
      core == 1 ~ "Treated_Core",
      b5   == 1 ~ "Treated_Band_5km",
      b10  == 1 ~ "Treated_Band_10km",
      b15  == 1 ~ "Treated_Band_15km",
      b20  == 1 ~ "Treated_Band_20km",
      TRUE      ~ "Control"
    )
  )

# ------------- 밴드별 event-study DiD (TWFE, ref = -1 = 2012) ---------------- #
m_es <- feols(
  outcome ~
    i(rel_time, core, ref = -1) +
    i(rel_time, b5,   ref = -1) +
    i(rel_time, b10,  ref = -1) +
    i(rel_time, b15,  ref = -1) +
    i(rel_time, b20,  ref = -1)
  | unit_id + year,
  data = did_df, vcov = ~ unit_id
)
etable(m_es)

# --------------------------- 계수 정리 -------------------------------------- #
es_df <- broom::tidy(m_es, conf.int = TRUE) %>%
  filter(grepl("^rel_time::", term)) %>%
  mutate(
    rel_time = as.integer(sub("^rel_time::(-?[0-9]+):.*$", "\\1", term)),
    band = case_when(
      grepl(":core$", term) ~ "Treated_Core",
      grepl(":b5$",   term) ~ "Treated_Band_5km",
      grepl(":b10$",  term) ~ "Treated_Band_10km",
      grepl(":b15$",  term) ~ "Treated_Band_15km",
      grepl(":b20$",  term) ~ "Treated_Band_20km",
      TRUE ~ NA_character_
    ),
    year = TREAT_YEAR + rel_time
  ) %>%
  filter(!is.na(band)) %>%
  mutate(band = factor(band, levels = treated_cols))

# ------------------------------- 그림 --------------------------------------- #
COL_DID <- "#0072B2"
p_did_es <- es_df %>%
  ggplot(aes(x = year, y = estimate)) +
  theme(aspect_ratio = 0.2) +
  scale_x_continuous(breaks = c(1997, 2007, 2013)) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = COL_DID, alpha = 0.25) +
  geom_line(color = COL_DID, linewidth = 1.1) +
  geom_vline(xintercept = 2013, linewidth = 0.7) +
  facet_wrap(~ band, nrow = 1, scales = "free_y") +
  labs(x = "Years", y = "DiD Event-study (ref = -1)") +
  theme_bw(base_size = 14) + theme(panel.grid.minor = element_blank())

ggsave(outf("did_event_study.png"), p_did_es, width = 14, height = 3.6, dpi = 500)
message("DiD 완료: ", outf("did_event_study.png"))
