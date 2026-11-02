# ==============================================================================
# SCRIPT ULTIMATE
# ==============================================================================

rm(list = ls())
set.seed(123)

# 1. LOAD LIBRARIES
library(wavethresh); library(breakfast); library(changepoint)
library(cpm); library(waveslim); library(dplyr)
library(ggplot2); library(tidyr); library(writexl)

# Pastikan select()/filter() selalu memakai versi dplyr, meskipun ada paket lain
select <- dplyr::select
filter <- dplyr::filter

# --- CONFIG ---
n <- 1024
n_iter <- 1000
data_types <- c("Big Jump", "Close Jumps", "Weak Signal", "Asymmetric",
                "Short Burst", "Trend + Jump", "Noise Variation", "Complex Mix")

sigma_fixed <- 0.5                     # sigma tetap untuk skenario 1-6
noise_variation_sigmas <- c(1.0, 2.0)  # khusus skenario "Noise Variation" 

# --- FUNGSI EVALUASI (NO DOUBLE COUNTING) ---
eval_met <- function(true_cp, det_cp, tol = 3) {
  tp <- 0
  if (length(det_cp) == 0) return(c(Sensitivity = 0, F1 = 0, Eta_Error = -length(true_cp)))

  matched <- rep(FALSE, length(det_cp))
  for (t in true_cp) {
    dist <- abs(det_cp - t)
    idx <- which(dist <= tol & !matched)
    if (length(idx) > 0) {
      tp <- tp + 1
      matched[idx[1]] <- TRUE
    }
  }
  rec <- tp / length(true_cp)
  prec <- tp / length(det_cp)
  f1 <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  return(c(Sensitivity = rec, F1 = f1, Eta_Error = length(det_cp) - length(true_cp)))
}

# --- FUNGSI REKONSTRUKSI SINYAL (untuk metrik MSE) ---
# Membangun estimasi sinyal piecewise-constant (y_hat) dari titik perubahan
# yang terdeteksi, dengan mean tiap segmen dihitung dari data observasi (y).
est_signal <- function(y, cp, n) {
  cp <- sort(unique(cp))
  cp <- cp[cp > 0 & cp < n]   # buang titik di luar batas wajar (termasuk cpts() PELT yang menyertakan n)
  bounds <- c(0, cp, n)
  yhat <- numeric(n)
  for (k in 1:(length(bounds) - 1)) {
    idx <- (bounds[k] + 1):bounds[k + 1]
    yhat[idx] <- mean(y[idx])
  }
  yhat
}

results_list <- list(); fingerprint_list <- list()
true_signal_list <- list(); obs_data_list <- list(); approx_signal_list <- list()
methods_names <- c("Haar_Basic", "Haar_FDR", "PELT", "MODWT", "TGUH", "WBS", "CUSUM")
counter <- 1; fp_counter <- 1; sig_counter <- 1; approx_counter <- 1

# --- 2. LOOP SIMULASI ---
for (type in data_types) {
  sigmas_to_run <- if (type == "Noise Variation") noise_variation_sigmas else sigma_fixed
  for (sig in sigmas_to_run) {
    cat("\n>>> Memproses:", type, "| Sigma:", sig, "\n")
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
    hit_counts <- matrix(0, nrow = n, ncol = 7)
    colnames(hit_counts) <- methods_names

    for (i in 1:n_iter) {
      setTxtProgressBar(pb, i)

      # Generator Sinyal
      if (type %in% c("Big Jump", "Noise Variation")) { f <- c(rep(0, 300), rep(4, 724)); true_cp <- c(301)
      } else if (type == "Close Jumps") { f <- c(rep(0, 400), rep(3, 20), rep(0, 604)); true_cp <- c(401, 421)
      } else if (type == "Weak Signal") { f <- c(rep(0, 500), rep(0.8, 524)); true_cp <- c(501)
      } else if (type == "Asymmetric") { f <- c(rep(0, 250), rep(3, 300), rep(2, 474)); true_cp <- c(251, 551)
      } else if (type == "Short Burst") { f <- c(rep(0, 495), rep(5, 10), rep(0, 519)); true_cp <- c(496, 506)
      } else if (type == "Trend + Jump") { f <- c(rep(0, 300), rep(2, 400), rep(0, 324)) + (0.005 * (1:n)); true_cp <- c(301, 701)
      } else if (type == "Complex Mix") {
        # Gabungan pola naik-turun tidak beraturan, menguji ketangguhan algoritma
        # pada rangkaian titik perubahan berdekatan dengan arah & amplitudo
        # yang saling berlawanan -- lebih menyerupai kondisi data riil.
        f <- c(rep(0,   100),   # baseline
               rep(3,   150),   # naik (big jump)                cp: 101
               rep(6,    15),   # naik lagi, short spike         cp: 251
               rep(1,    60),   # turun drastis di bawah awal    cp: 266
               rep(4,   200),   # naik lagi                      cp: 326
               rep(8,    10),   # short burst ekstrem             cp: 526
               rep(4.5,  80),   # turun sebagian                 cp: 536
               rep(2,    30),   # turun lagi                     cp: 616
               rep(5,   120),   # naik                            cp: 646
               rep(5.3,  25),   # weak shift kecil                cp: 766
               rep(1.5,  90),   # turun drastis                   cp: 791
               rep(3,   144) + 0.006 * (1:144))  # naik + tren lokal  cp: 881
        true_cp <- c(101, 251, 266, 326, 526, 536, 616, 646, 766, 791, 881)
      }

      err <- rnorm(n, 0, sig)
      y <- f + err

      if (i == 1) {
        true_signal_list[[sig_counter]] <- data.frame(Skenario = type, Sigma = sig, Index = 1:n, Nilai = f, Tipe = "True Signal")
        obs_data_list[[sig_counter]] <- data.frame(Skenario = type, Sigma = sig, Index = 1:n, Nilai = y)
        sig_counter <- sig_counter + 1
      }

      # Deteksi
      wd_y <- wd(y, filter.number = 1, family = "DaubExPhase")
      cp_h <- which(diff(wr(threshold(wd_y, policy = "universal"))) != 0) + 1
      cp_f <- suppressWarnings(which(diff(wr(threshold(wd_y, policy = "fdr", alpha = 0.05))) != 0) + 1)
      cp_p <- cpts(cpt.mean(y, method = "PELT"))

      # MODWT (MAD Threshold)
      W <- modwt(y, wf = "haar", n.levels = 1)
      sigma_mad <- median(abs(W[[1]])) / 0.6745
      thr_univ <- sigma_mad * sqrt(2 * log(n))
      cp_m <- which(abs(W[[1]]) > thr_univ)

      cp_t <- model.ic(sol.tguh(y))$cpts
      cp_w <- model.ic(sol.wbs(y))$cpts

      # CUSUM
      cp_c <- processStream(y, cpmType = "Student")$changePoints

      m_list <- list(Haar_Basic = cp_h, Haar_FDR = cp_f, PELT = cp_p, MODWT = cp_m, TGUH = cp_t, WBS = cp_w, CUSUM = cp_c)

      for (m in methods_names) {
        cp_est <- m_list[[m]]
        res <- eval_met(true_cp, cp_est)
        yhat_m <- est_signal(y, cp_est, n)
        mse_m <- mean((yhat_m - f)^2)

        results_list[[counter]] <- data.frame(Skenario = type, Sigma = sig, Metode = m,
                                               Sensitivity = res[1], F1 = res[2],
                                               Eta_Error = res[3], MSE = mse_m)
        counter <- counter + 1

        det <- cp_est; det <- det[det > 0 & det <= n]
        hit_counts[det, m] <- hit_counts[det, m] + 1

        if (i == 1) {
          approx_signal_list[[approx_counter]] <- data.frame(Skenario = type, Sigma = sig,
                                                               Index = 1:n, Nilai = yhat_m, Tipe = m)
          approx_counter <- approx_counter + 1
        }
      }
    }
    close(pb)
    fp_temp <- as.data.frame(hit_counts / n_iter); fp_temp$Index <- 1:n; fp_temp$Skenario <- type; fp_temp$Sigma <- sig
    fingerprint_list[[fp_counter]] <- fp_temp; fp_counter <- fp_counter + 1
  }
}

# --- 3. OUTPUT CONSOLE & TABEL (EXCEL) ---
results_raw <- bind_rows(results_list)
final_summary <- results_raw %>% group_by(Skenario, Sigma, Metode) %>%
  summarise(Mean_Sens = round(mean(Sensitivity), 4), Mean_F1 = round(mean(F1), 4),
            Eta_Under = sum(Eta_Error < 0), Eta_Exact = sum(Eta_Error == 0), Eta_Over = sum(Eta_Error > 0),
            E_Eta_Error = round(mean(Eta_Error), 4), Var_Eta = round(var(Eta_Error), 4),
            Mean_MSE = round(mean(MSE), 6), .groups = 'drop')

cat("\n\n=== HASIL REKAPITULASI FINAL (CUSUM: processStream) ===\n")
print(as.data.frame(final_summary))

cat("\n>>> Menyimpan tabel ke 'HASIL_FINAL_TABEL_processStream.xlsx'...\n")
write_xlsx(final_summary, "HASIL_FINAL_TABEL_processStream.xlsx")

# --- WARNA (dipakai konsisten di semua PNG per metode) ---
my_colors <- c("True Signal" = "#333333", "Haar_Basic" = "#E41A1C", "Haar_FDR" = "#377EB8",
               "PELT" = "#4DAF4A", "MODWT" = "#984EA3", "TGUH" = "#FF7F00", "WBS" = "#A65628", "CUSUM" = "#F781BF")

sanitize_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# --- 3.5 VISUALISASI BARU: DATA ASLI (SINYAL ASLI + NOISE) PER SKENARIO ---
# 1 PNG per kombinasi skenario+sigma, menampilkan garis sinyal asli f(t) dan
# sebaran data observasi (y) dari realisasi noise pada iterasi pertama.
da_root <- "Plot_Data_Asli_processStream"
cat("\n>>> Menyimpan plot data asli (sinyal + noise) ke folder '", da_root, "/'...\n", sep = "")
dir.create(da_root, showWarnings = FALSE)

true_long <- bind_rows(true_signal_list)   # Skenario, Sigma, Index, Nilai, Tipe = "True Signal"
obs_long <- bind_rows(obs_data_list)       # Skenario, Sigma, Index, Nilai (data y observasi)

da_combos <- true_long %>% select(Skenario, Sigma) %>% distinct()

for (i in 1:nrow(da_combos)) {
  scen <- da_combos$Skenario[i]; sig <- da_combos$Sigma[i]

  df_true_i <- true_long %>% filter(Skenario == scen, Sigma == sig)
  df_obs_i <- obs_long %>% filter(Skenario == scen, Sigma == sig)

  p_da <- ggplot() +
    geom_line(data = df_obs_i, aes(x = Index, y = Nilai), color = "#AAB8D4", size = 0.3, alpha = 0.8) +
    geom_step(data = df_true_i, aes(x = Index, y = Nilai), color = "#1B2A6B", size = 0.8) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank()) +
    labs(title = paste0("DATA ASLI: ", scen, " (Sigma = ", sig, ")"), x = "t", y = "Nilai")

  fname <- paste0(sanitize_name(scen), "_sigma_", sig, ".png")
  ggsave(file.path(da_root, fname), plot = p_da, width = 9, height = 2.5, dpi = 150, bg = "white")
}

# --- 4. VISUALISASI: FINGERPRINT (PROPORSI DETEKSI) -- 1 PNG BERSIH PER METODE PER SKENARIO ---
fp_root <- "Plot_Fingerprint_processStream"
cat("\n>>> Menyimpan plot fingerprint (PNG per metode) ke folder '", fp_root, "/'...\n", sep = "")
dir.create(fp_root, showWarnings = FALSE)

prop_long <- bind_rows(fingerprint_list) %>%
  pivot_longer(cols = all_of(methods_names), names_to = "Tipe", values_to = "Nilai")

fp_combos <- prop_long %>% select(Skenario, Sigma) %>% distinct()

for (i in 1:nrow(fp_combos)) {
  scen <- fp_combos$Skenario[i]; sig <- fp_combos$Sigma[i]
  sub_folder <- file.path(fp_root, paste0(sanitize_name(scen), "_sigma_", sig))
  dir.create(sub_folder, showWarnings = FALSE, recursive = TRUE)

  for (m in methods_names) {
    df_m <- prop_long %>% filter(Skenario == scen, Sigma == sig, Tipe == m)
    df_m$Tipe <- factor(df_m$Tipe, levels = m)

    p <- ggplot(df_m, aes(x = Index, y = Nilai)) +
      geom_line(color = my_colors[[m]], size = 0.5) +
      facet_grid(Tipe ~ .) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      theme_bw(base_size = 9) +
      theme(strip.text.y = element_text(angle = 0, face = "bold", size = 9),
            strip.background = element_rect(fill = "white"),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            plot.margin = margin(3, 3, 3, 3))

    ggsave(file.path(sub_folder, paste0(m, ".png")), plot = p,
           width = 8, height = 1.3, dpi = 150, bg = "white")
  }
}

# --- 5. VISUALISASI: HAMPIRAN SINYAL -- 1 PNG BERSIH PER METODE PER SKENARIO (+ sinyal asli) ---
hs_root <- "Plot_Hampiran_Sinyal_processStream"
cat("\n>>> Menyimpan plot hampiran sinyal (PNG per metode) ke folder '", hs_root, "/'...\n", sep = "")
dir.create(hs_root, showWarnings = FALSE)

approx_long <- bind_rows(approx_signal_list)   # Skenario, Sigma, Index, Nilai, Tipe (nama metode)
# (true_long dan obs_long sudah dihitung di Section 3.5)

hs_combos <- true_long %>% select(Skenario, Sigma) %>% distinct()

for (i in 1:nrow(hs_combos)) {
  scen <- hs_combos$Skenario[i]; sig <- hs_combos$Sigma[i]
  sub_folder <- file.path(hs_root, paste0(sanitize_name(scen), "_sigma_", sig))
  dir.create(sub_folder, showWarnings = FALSE, recursive = TRUE)

  df_obs <- obs_long %>% filter(Skenario == scen, Sigma == sig)

  # -- Panel referensi: sinyal asli (True Signal) --
  df_true <- true_long %>% filter(Skenario == scen, Sigma == sig)
  df_true$Tipe <- factor(df_true$Tipe, levels = "True Signal")

  p_true <- ggplot() +
    geom_point(data = df_obs, aes(x = Index, y = Nilai), color = "grey70", size = 0.3, alpha = 0.5) +
    geom_step(data = df_true, aes(x = Index, y = Nilai), color = my_colors[["True Signal"]], size = 0.7) +
    facet_grid(Tipe ~ .) +
    theme_bw(base_size = 9) +
    theme(strip.text.y = element_text(angle = 0, face = "bold", size = 9),
          strip.background = element_rect(fill = "white"),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          plot.margin = margin(3, 3, 3, 3))

  ggsave(file.path(sub_folder, "True_Signal.png"), plot = p_true,
         width = 8, height = 1.3, dpi = 150, bg = "white")

  # -- Panel tiap metode (dengan garis putus-putus sinyal asli sebagai referensi) --
  df_true_ref <- df_true %>% select(Index, Nilai)   # tanpa kolom Tipe, supaya tidak bikin panel baru

  for (m in methods_names) {
    df_m <- approx_long %>% filter(Skenario == scen, Sigma == sig, Tipe == m)
    df_m$Tipe <- factor(df_m$Tipe, levels = m)

    p_m <- ggplot() +
      geom_step(data = df_true_ref, aes(x = Index, y = Nilai), color = "grey35", linetype = "dashed", size = 0.5) +
      geom_step(data = df_m, aes(x = Index, y = Nilai), color = my_colors[[m]], size = 0.7) +
      facet_grid(Tipe ~ .) +
      theme_bw(base_size = 9) +
      theme(strip.text.y = element_text(angle = 0, face = "bold", size = 9),
            strip.background = element_rect(fill = "white"),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            plot.margin = margin(3, 3, 3, 3))

    ggsave(file.path(sub_folder, paste0(m, ".png")), plot = p_m,
           width = 8, height = 1.3, dpi = 150, bg = "white")
  }
}

cat("\n>>> SELESAI! Cek folder kerjamu:\n")
cat("    -", "HASIL_FINAL_TABEL_processStream.xlsx", "(tabel metrik)\n")
cat("    -", fp_root, "/<Skenario>_sigma_<sigma>/<Metode>.png (fingerprint per metode)\n")
cat("    -", hs_root, "/<Skenario>_sigma_<sigma>/<Metode>.png + True_Signal.png (hampiran sinyal)\n")
