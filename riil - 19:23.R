# ==============================================================================
# SCRIPT GLOBAL: DETEKSI JAKARTA (2019 - 2023)
# 4 PANEL (Data Asli, PELT, WBS, TGUH)
# ==============================================================================

rm(list = ls())
library(dplyr)
library(lubridate)
library(imputeTS)
library(ggplot2)
library(readr)
library(changepoint) 
library(wbs)         
library(breakfast)   

# 1. PERSIAPAN DATA (Fokus 2019 - 2023)
df_raw <- read_csv("jakarta-air-quality.csv")
colnames(df_raw) <- trimws(colnames(df_raw))
df_raw$date <- ymd(trimws(df_raw$date))
df_raw$pm25 <- as.numeric(trimws(df_raw$pm25))

kalender <- data.frame(date = seq(as.Date("2019-01-01"), as.Date("2023-12-31"), by = "day"))
df_clean <- left_join(kalender, df_raw[, c("date", "pm25")], by = "date")

impute_stochastic_local <- function(y, radius = 100) {
  idx_na <- which(is.na(y))
  if(length(idx_na) == 0) return(y)
  interp_base <- na_interpolation(y, option = "linear")
  n_data <- length(y)
  local_sd <- numeric(length(idx_na))
  for (i in seq_along(idx_na)) {
    start_idx <- max(1, idx_na[i] - radius)
    end_idx <- min(n_data, idx_na[i] + radius)
    sd_val <- sd(y[start_idx:end_idx], na.rm = TRUE)
    if (is.na(sd_val)) sd_val <- sd(y, na.rm = TRUE) * 0.5
    local_sd[i] <- sd_val
  }
  set.seed(10122054) 
  y_stochastic <- interp_base
  y_stochastic[idx_na] <- interp_base[idx_na] + rnorm(length(idx_na), 0, local_sd)
  y_stochastic[y_stochastic < 5] <- 5 
  return(y_stochastic)
}

y_clean <- impute_stochastic_local(df_clean$pm25, radius = 100)
df_clean$pm25_stochastic <- y_clean

# 2. EKSEKUSI 3 ALGORITMA UTAMA
cp_pelt <- cpts(cpt.mean(y_clean, method = "PELT", penalty = "BIC"))
cp_wbs <- changepoints(wbs(y_clean), penalty = "ssic.penalty")$cpt.ic$ssic.penalty
cp_tguh <- model.ic(sol.tguh(y_clean))$cpts 

# 3. PENYUSUNAN PLOT
level_urutan <- c("Data Asli", "PELT", "WBS", "TGUH")
df_cp <- data.frame(
  Metode = factor(c(rep("PELT", length(cp_pelt)), rep("WBS", length(cp_wbs)), rep("TGUH", length(cp_tguh))), levels = level_urutan),
  Tanggal_CP = df_clean$date[c(cp_pelt, cp_wbs, cp_tguh)]
)

df_plot <- bind_rows(
  df_clean %>% mutate(Metode = "Data Asli"),
  df_clean %>% mutate(Metode = "PELT"),
  df_clean %>% mutate(Metode = "WBS"),
  df_clean %>% mutate(Metode = "TGUH")
) %>% mutate(Metode = factor(Metode, levels = level_urutan))

ggplot(df_plot, aes(x = date, y = pm25_stochastic)) +
  geom_line(color = "gray30", linewidth = 0.3) +
  geom_vline(data = subset(df_cp, !is.na(Tanggal_CP)), 
             aes(xintercept = as.numeric(Tanggal_CP), color = Metode), linewidth = 0.6, alpha = 0.8) +
  scale_color_manual(values = c("PELT" = "#4DAF4A", "WBS" = "#E41A1C", "TGUH" = "#FF7F00")) +
  facet_grid(Metode ~ .) +
  labs(title = "Deteksi Titik Perubahan PM 2.5 Jakarta (2019-2023)",
       subtitle = "Komparasi Kinerja Algoritma PELT, WBS, dan TGUH",
       x = "Tahun", y = "Konsentrasi PM 2.5 (\u03bcg/m\u00b3)") +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 10))