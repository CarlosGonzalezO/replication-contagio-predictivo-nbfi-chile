# =============================================================================
# Script: 03_figuras_dy_BCCh.R
# Autor:  Carlos Gonzalez
# 
# Genera todas las figuras de la tesis a partir del bootstrap pre-computado.
# Requiere haber corrido previamente 02_bootstrap_dy_BCCh.R que genera el RDS.
# No obstante para efectos de replicación en carpeta se halla el RDS.
# =============================================================================
# =============================================================================
# PASO 0 - RUTA
# =============================================================================
Sys.setlocale("LC_ALL", "es_ES.UTF-8")
rm(list = ls())

# Por favor, setear ruta en "ruta" y luego ejecutar código.
library(fs)
ruta <- path("../..")
ruta_input  <- ruta / "1_datos/1_clean_data"
ruta_output <- ruta / "3.2_resultado_dy/salidas_bootstrap"

# =============================================================================
# CARGAR BOOTSTRAP PRE-COMPUTADO
# =============================================================================
library(tidyverse)
library(lubridate)
cat("Cargando bootstrap pre-computado...\n")
boot <- readRDS(file.path(ruta_output, "bcch_bootstrap_rolling.rds"))

# Asignar a variables locales (más cómodo para escribir el código)
boot_tci_rw <- boot$boot_tci_rw
boot_net_rw <- boot$boot_net_rw
boot_C_rw   <- boot$boot_C_rw
df_net_boot <- boot$df_net_boot
df_tci_boot <- boot$df_tci_boot
TCI_obs_rw  <- boot$TCI_obs_rw
NET_obs_rw  <- boot$NET_obs_rw
CT_rw_obs   <- boot$CT_rw_obs
fechas_rw   <- boot$fechas_rw
nombre_vars <- boot$nombre_vars
B_rolling   <- boot$B_rolling
WINDOW      <- boot$WINDOW
K           <- boot$K
ALPHA_90    <- boot$ALPHA_90
ALPHA_68    <- boot$ALPHA_68

n_dates <- length(fechas_rw)
q90_lo  <- ALPHA_90 / 2; q90_hi <- 1 - ALPHA_90 / 2
q68_lo  <- ALPHA_68 / 2; q68_hi <- 1 - ALPHA_68 / 2

cat(sprintf("Cargado: B=%d, n_dates=%d, K=%d, W=%d\n", 
            B_rolling, n_dates, K, WINDOW))

tema_g <- theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "gray95", color = NA),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 40, hjust = 1, color = "gray30"),
        axis.text.y = element_text(color = "gray30"),
        panel.grid.minor = element_blank(),
        #plot.title = element_text(face = "bold", size = 12),
        #plot.subtitle = element_text(color = "gray40", size = 9)
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(family = "Times New Roman", face = "bold", size = 12),
        plot.subtitle = element_text(family = "Times New Roman", color = "gray30", size = 9),)

#plot.title = element_text(face = "bold", size = 12),
theme(panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray92", linewidth = 0.3))
#plot.subtitle = element_text(color = "gray40", size = 9)
theme(panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.4))

# Orden de narración
orden_narrativo <- c("FP", "Bancos", "BancoCentral", "Seguros", 
                     "FMM", "FNM", "OFIs", "RestoMundo")
df_net_boot <- df_net_boot %>% 
  mutate(sector = factor(sector, levels = orden_narrativo))

eventos <- function() list(
  #geom_vline(xintercept = as.Date("2008-09-01"), color = "gray40", linetype = "dotted"),
  geom_vline(xintercept = as.Date("2019-10-01"), color = "gray40", linetype = "dotted"),
  geom_vline(xintercept = as.Date(c("2020-07-01","2020-12-01","2021-04-01")),
             color = "red", linetype = "dotted", linewidth = 0.4),
  geom_vline(xintercept = as.Date("2020-07-01"), 
             color = "red", linetype = "dotted", linewidth = 0.5),
  # annotate("text", x = as.Date("2020-07-01"), y = -Inf, 
  #          label = "Inicio Retiros", vjust = -0.5, hjust = -0.1,
  #          size = 2.8, color = "gray20", family = "Times New Roman"),
  annotate("rect", xmin = as.Date("2020-07-01"), xmax = as.Date("2021-12-31"),
           ymin = -Inf, ymax = Inf, alpha = 0.10, fill = "red"))

# TCI: dual band + bias-corrected
p_tci <- ggplot(df_tci_boot, aes(x = fecha)) +
  geom_ribbon(aes(ymin = tci_lo90, ymax = tci_hi90), fill = "steelblue", alpha = 0.18) +
  geom_ribbon(aes(ymin = tci_lo68, ymax = tci_hi68), fill = "steelblue", alpha = 0.35) +
  #geom_line(aes(y = tci_obs), color = "gray50", linewidth = 0.5, linetype = "dashed") +
  geom_line(aes(y = tci_BC),  color = "steelblue4", linewidth = 0.9) +
  geom_hline(yintercept = 36.87, color = "gray30", linetype = "dashed", linewidth = 0.4) +
  annotate("text", x = as.Date("2013-01-01"), y = 38.5,
           label = "Baseline Pre-Estallido (36.9%)",
           color = "gray30", size = 3, hjust = 0, family = "Times New Roman") +
  eventos() +
  annotate("text", x = as.Date("2021-04-01"), y = Inf, label = "Periodo\nRetiros",
           vjust = 1.5, hjust = 0.5, size = 3, color = "darkred", fontface = "bold", family = "Times New Roman") +
  scale_x_date(
    date_breaks = "2 years",
    date_minor_breaks = "1 year",
    labels = function(x) ifelse(is.na(x), "",
                                paste0(lubridate::year(x))),  # solo año
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = "gray20")) +
  labs(title = "Índice de Conectividad Total (TCI)",
       subtitle = "Conectividad financiera, 2011T3-2025T3. Bandas al 68% y 90%.",
       x = "Periodo", y = "TCI (%)") +
  tema_g
ggsave(file.path(ruta_output, "figs", "bcch_15_tci_bootstrap_BC.png"),
       p_tci, width = 10, height = 6, dpi = 150, bg = "white")

# NET: dual band + bias-corrected
p_net <- ggplot(df_net_boot, aes(x = fecha)) +
  geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "steelblue", alpha = 0.15) +
  geom_ribbon(aes(ymin = lo68, ymax = hi68), fill = "steelblue", alpha = 0.32) +
  #geom_line(aes(y = obs), color = "gray55", linewidth = 0.4, linetype = "dashed") +
  geom_line(aes(y = BC), color = "steelblue4", linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  facet_wrap(~ sector, scales = "free_y", ncol = 2) +
  eventos() +
  annotate("text", x = as.Date("2021-04-01"), y = Inf, label = "Periodo\nRetiros",
           vjust = 1.5, hjust = 0.5, size = 3, color = "darkred", fontface = "bold", family = "Times New Roman") +
  scale_x_date(
    date_breaks = "2 years",
    date_minor_breaks = "1 year",
    labels = function(x) ifelse(is.na(x), "",
                                paste0(lubridate::year(x))),  # solo año
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = "gray20")) +
  labs(title = "Connectividad neta direccional por sector, 2011T3-2025T3",
       subtitle = "Valores positivos indican transmisor neto de shocks.",
       x = "Periodo", y = "Connectividad neta") +
  tema_g
ggsave(file.path(ruta_output, "figs", "bcch_16_net_bootstrap_BC.png"),
       p_net, width = 10, height = 9, dpi = 150, bg = "white")

# =============================================================================
# (4) REGIMENES Y CONTRASTES (con IC basico)
# =============================================================================
cat("\n--- (4) REGIMENES Y CONTRASTES ---\n")

regimenes <- tibble::tribble(
  ~regimen,         ~inicio,         ~fin,
  "Pre-Estallido",  "2012-01-01",    "2019-09-30",   # ~33 obs
  "Estallido",      "2019-10-01",    "2020-03-30",   # ~3  obs
  "Retiros",        "2020-04-01",    "2021-06-30",   # ~6  obs (considero trimestre inmediatamente anterior y posterior)
  "Post-Retiros",   "2022-01-01",    "2025-09-30"    # ~15 obs
)
#promedia TCI/NET de esos índices, tanto en la data observada como en cada réplica bootstrap.
#La identificación de regímenes se basa en un dating exógeno de eventos políticos y de política pública chilenos. El régimen Pre-Estallido (2011T3–2019T3) corresponde al período de normalidad institucional posterior a la recuperación de la crisis financiera global. El régimen Estallido (2019T4–2020T2) inicia con el comienzo del estallido social en octubre de 2019 y se extiende hasta justo antes del primer retiro de fondos previsionales, capturando además los meses iniciales de la pandemia. El régimen Retiros (2020T3–2021T4) abarca el período durante el cual se materializaron los tres retiros excepcionales de fondos previsionales (julio 2020, diciembre 2020 y abril 2021), incluyendo los trimestres inmediatamente posteriores para capturar persistencia de corto plazo. El régimen Post-Retiros (2022T1–2025T3) corresponde al período de normalización posterior. La definición se realiza exógenamente al comportamiento del TCI, evitando data snooping en la selección de fechas.



get_idx <- function(reg) {
  r <- regimenes[regimenes$regimen == reg, ]
  which(fechas_rw >= r$inicio & fechas_rw <= r$fin)
}
for (rg in regimenes$regimen) cat(sprintf("  %s: %d obs\n", rg, length(get_idx(rg))))

# Medias por regimen, observado y bootstrap
tci_obs_reg <- sapply(regimenes$regimen, function(r) {
  idx <- get_idx(r); if (length(idx) > 0) mean(TCI_obs_rw[idx], na.rm = TRUE) else NA_real_
})
tci_boot_reg <- sapply(regimenes$regimen, function(r) {
  idx <- get_idx(r)
  if (length(idx) > 0) colMeans(boot_tci_rw[idx, , drop = FALSE], na.rm = TRUE)
  else rep(NA_real_, B_rolling)
})
colnames(tci_boot_reg) <- regimenes$regimen

net_obs_reg  <- matrix(NA_real_, K, nrow(regimenes), dimnames = list(nombre_vars, regimenes$regimen))
net_boot_reg <- array(NA_real_, c(B_rolling, K, nrow(regimenes)),
                      dimnames = list(NULL, nombre_vars, regimenes$regimen))
for (r in 1:nrow(regimenes)) {
  idx <- get_idx(regimenes$regimen[r]); if (length(idx) == 0) next
  for (s in 1:K) {
    net_obs_reg[s, r] <- mean(NET_obs_rw[idx, s], na.rm = TRUE)
    net_boot_reg[, s, r] <- apply(boot_net_rw[idx, s, , drop = FALSE], 3, mean, na.rm = TRUE)
  }
}

#funcion qfun
qfun <- function(x, p) quantile(x, probs = p, na.rm = TRUE, names = FALSE)

#funcion ic_basico_escalar
ic_basico_escalar <- function(theta_hat, boot_vec) {
  theta_bar <- mean(boot_vec, na.rm = TRUE)
  bias      <- theta_bar - theta_hat
  list(
    theta_bar = theta_bar,
    bias      = bias,
    theta_BC  = theta_hat - bias,
    lo90 = 2*theta_hat - qfun(boot_vec, q90_hi),
    hi90 = 2*theta_hat - qfun(boot_vec, q90_lo),
    lo68 = 2*theta_hat - qfun(boot_vec, q68_hi),
    hi68 = 2*theta_hat - qfun(boot_vec, q68_lo)
  )
}

# IC basico para medias por regimen
df_regime_tci <- bind_rows(lapply(regimenes$regimen, function(r) {
  ic <- ic_basico_escalar(tci_obs_reg[r], tci_boot_reg[, r])
  tibble(regimen = r, obs = tci_obs_reg[r],
         BC = ic$theta_BC,
         ic68_lo = ic$lo68, ic68_hi = ic$hi68,
         ic90_lo = ic$lo90, ic90_hi = ic$hi90)
})) %>% mutate(regimen = factor(regimen, levels = regimenes$regimen))
write_csv(df_regime_tci, file.path(ruta_output, "bcch_18_regime_means_tci.csv"))

df_regime_net <- bind_rows(lapply(seq_len(nrow(regimenes)), function(r_i) {
  bind_rows(lapply(seq_len(K), function(s) {
    obs <- net_obs_reg[s, r_i]
    ic  <- ic_basico_escalar(obs, net_boot_reg[, s, r_i])
    tibble(regimen = regimenes$regimen[r_i], sector = nombre_vars[s],
           obs = obs, BC = ic$theta_BC,
           ic68_lo = ic$lo68, ic68_hi = ic$hi68,
           ic90_lo = ic$lo90, ic90_hi = ic$hi90)
  }))
})) %>% mutate(sector = factor(sector, levels = nombre_vars),
               regimen = factor(regimen, levels = regimenes$regimen))
write_csv(df_regime_net, file.path(ruta_output, "bcch_19_regime_means_net.csv"))

# Contrastes (diferencia A - B). Para diferencias, el sesgo se cancela
# parcialmente; el IC basico aplica igualmente.
orden_contrastes_tci <- c("Estallido - Pre-Estallido",
                          "Retiros - Pre-Estallido",
                          "Post-Retiros - Pre-Estallido",
                          "Retiros - Estallido",
                          "Post-Retiros - Retiros")


hacer_contraste <- function(boot_A, boot_B, obs_A, obs_B) {
  ic <- ic_basico_escalar(obs_A - obs_B, boot_A - boot_B)
  diff_dist <- boot_A - boot_B
  p_val <- 2 * min(mean(diff_dist <= 0, na.rm = TRUE),
                   mean(diff_dist >= 0, na.rm = TRUE))
  list(obs = obs_A - obs_B, BC = ic$theta_BC,
       ic68_lo = ic$lo68, ic68_hi = ic$hi68,
       ic90_lo = ic$lo90, ic90_hi = ic$hi90,
       p_value = p_val)
}

contrastes <- tibble::tribble(
  ~A,              ~B,
  "Estallido",     "Pre-Estallido",
  "Retiros",       "Pre-Estallido",
  "Retiros",       "Estallido",
  "Post-Retiros",  "Pre-Estallido",
  "Post-Retiros",  "Retiros"
)

contrast_tci <- bind_rows(lapply(1:nrow(contrastes), function(k) {
  A <- contrastes$A[k]; B_ <- contrastes$B[k]
  ct <- hacer_contraste(tci_boot_reg[, A], tci_boot_reg[, B_],
                        tci_obs_reg[A], tci_obs_reg[B_])
  tibble(contraste = paste0(A, " - ", B_), metric = "TCI", sector = NA_character_,
         obs = ct$obs, BC = ct$BC,
         ic68_lo = ct$ic68_lo, ic68_hi = ct$ic68_hi,
         ic90_lo = ct$ic90_lo, ic90_hi = ct$ic90_hi,
         p_value = ct$p_value,
         sig_68 = !(ct$ic68_lo <= 0 & ct$ic68_hi >= 0),
         sig_90 = !(ct$ic90_lo <= 0 & ct$ic90_hi >= 0))
}))

contrast_tci <- contrast_tci %>% 
  mutate(contraste = factor(contraste, levels = orden_contrastes_tci))

contrast_tci <- contrast_tci %>%
  mutate(contraste = str_replace(contraste, " - ", " vs "))

contrast_net <- bind_rows(lapply(1:nrow(contrastes), function(k) {
  A <- contrastes$A[k]; B_ <- contrastes$B[k]
  bind_rows(lapply(1:K, function(s) {
    sname <- nombre_vars[s]
    ct <- hacer_contraste(net_boot_reg[, s, A], net_boot_reg[, s, B_],
                          net_obs_reg[s, A], net_obs_reg[s, B_])
    tibble(contraste = paste0(A, " - ", B_), metric = "NET", sector = sname,
           obs = ct$obs, BC = ct$BC,
           ic68_lo = ct$ic68_lo, ic68_hi = ct$ic68_hi,
           ic90_lo = ct$ic90_lo, ic90_hi = ct$ic90_hi,
           p_value = ct$p_value,
           sig_68 = !(ct$ic68_lo <= 0 & ct$ic68_hi >= 0),
           sig_90 = !(ct$ic90_lo <= 0 & ct$ic90_hi >= 0))
  }))
}))

df_contrastes <- bind_rows(contrast_tci, contrast_net) %>%
  mutate(stars = case_when(
    p_value < 0.01 ~ "***", p_value < 0.05 ~ "**",
    p_value < 0.10 ~ "*",   p_value < 0.32 ~ ".",   TRUE ~ ""
  ))
write_csv(df_contrastes, file.path(ruta_output, "bcch_20_regime_contrasts.csv"))

cat("\nContrastes TCI (con IC basico):\n")
print(contrast_tci %>%
        mutate(across(c(obs, BC, ic68_lo, ic68_hi, ic90_lo, ic90_hi), ~round(., 2)),
               p_value = round(p_value, 4)) %>%
        select(contraste, obs, BC, ic90_lo, ic90_hi, p_value, sig_90))

# --- Plot medias TCI por regimen ---
p_reg_tci <- ggplot(df_regime_tci, aes(x = regimen, y = BC)) +
  geom_col(fill = "steelblue", alpha = 0.4, width = 0.6) +
  geom_errorbar(aes(ymin = ic90_lo, ymax = ic90_hi), width = 0.15, linewidth = 0.4, color = "steelblue4") +
  geom_errorbar(aes(ymin = ic68_lo, ymax = ic68_hi), width = 0, linewidth = 2, color = "steelblue4") +
  geom_point(size = 2.5, color = "steelblue4") +
  #geom_text(aes(label = sprintf("%.1f", BC), y = ic90_hi), vjust = -0.5, size = 3.2, color = "steelblue4") +
  geom_text(aes(label = sprintf("%.1f", BC)), 
            hjust = -0.3, vjust = -0.6, size = 2.8, 
            color = "gray20", family = "Times New Roman") +
  labs(title = "TCI medio por regimen (bias-corrected) con IC basico 68% (grueso) y 90% (fino)",
       subtitle = sprintf("B=%d, W=%d", B_rolling, WINDOW), x = NULL, y = "TCI bias-corrected (%)") +
  tema_g
ggsave(file.path(ruta_output, "figs", "bcch_17_regime_means_tci_BC.png"),
       p_reg_tci, width = 9, height = 5.5, dpi = 150, bg = "white")

# --- Forest plots de contrastes ---
p_forest_tci <- ggplot(contrast_tci,
                       aes(x = BC, y = fct_rev(factor(contraste, levels = unique(contraste))))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = ic90_lo, xmax = ic90_hi, color = sig_90), height = 0, linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ic68_lo, xmax = ic68_hi, color = sig_90), height = 0, linewidth = 2) +
  geom_point(aes(color = sig_90), size = 2.8) +
  #geom_text(aes(label = sprintf("%.2f%s", BC,
  # case_when(p_value < 0.01 ~ "***", p_value < 0.05 ~ "**",
  #p_value < 0.10 ~ "*",   TRUE ~ ""))),
  #hjust = -0.3, vjust = -0.6, size = 3, color = "gray20") +
  geom_text(aes(label = sprintf("%.1f", BC)), 
            hjust = -0.3, vjust = -0.6, size = 2.8, 
            color = "gray20", family = "Times New Roman") +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "steelblue4"),
                     labels = c("FALSE" = "IC 90% contiene 0", "TRUE" = "Sig. 90%")) +
  labs(title = "Diferencias en TCI entre regímenes",
       subtitle = "Diferencia A − B en puntos porcentuales. Bandas: 68% (gruesa) y 90% (fina).",
       x = "", y = NULL, color = NULL) +
  tema_g
ggsave(file.path(ruta_output, "figs", "bcch_18_contrastes_tci_BC.png"),
       p_forest_tci, width = 9, height = 5, dpi = 150, bg = "white")

p_forest_net <- ggplot(contrast_net,
                       aes(x = BC, y = fct_rev(factor(sector, levels = nombre_vars)))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = ic90_lo, xmax = ic90_hi, color = sig_90), height = 0, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ic68_lo, xmax = ic68_hi, color = sig_90), height = 0, linewidth = 1.5) +
  geom_point(aes(color = sig_90), size = 2) +
  facet_wrap(~ factor(contraste, levels = unique(contraste)), 
             ncol = 2, scales = "free_x") +  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "steelblue4"),
                                                                labels = c("FALSE" = "IC 90% contiene 0", "TRUE" = "Sig. 90%")) +
  labs(title = "Contrastes en connectividad neta entre regímenes, por sector",
       subtitle = "Barra gruesa: IC basico 68%. Barra fina: IC 90%.",
       x = "Diferencia A − B en puntos porcentuales. Bandas: 68% (gruesa) y 90% (fina).", y = NULL, color = NULL) +
  tema_g
ggsave(file.path(ruta_output, "figs", "bcch_19_contrastes_net_BC.png"),
       p_forest_net, width = 9, height = 11, dpi = 150, bg = "white")

cat("\n===============================================================\n")
cat("BOOTSTRAP v5 COMPLETADO (con IC basico Hall 1992)\n")
cat("===============================================================\n")
cat("Linea solida en plots = estimado BIAS-CORRECTED.\n")
cat("Linea punteada gris   = estimado puntual sin corregir (para comparacion).\n")
cat("Banda gruesa interior = IC basico 68% (~ +/- 1 sigma).\n")
cat("Banda fina exterior   = IC basico 90%.\n")

# =============================================================================
# (5) FIGURAS PAIRWISE: FROM-hacia-FP y TO-desde-BancoCentral
# =============================================================================
cat("\n--- (5) PAIRWISE: FROM hacia FP y TO desde BancoCentral ---\n")

idx_FP   <- which(nombre_vars == "FP")
idx_BCCh <- which(nombre_vars == "BancoCentral")

# --- (a) FROM hacia FP: C[FP, j] para j != FP, por fecha ----------------------
# Observado: C_rw_obs[FP, j, t]
# Para FROM hacia FP: usar OBSERVADO + percentiles directos (no bias-correction)
from_FP_obs  <- t(sapply(1:n_dates, function(t) CT_rw_obs[idx_FP, , t]))
colnames(from_FP_obs) <- nombre_vars

# Percentiles directos sobre la distribucion bootstrap
from_FP_lo90B <- apply(boot_C_rw[, idx_FP, , ], c(1, 2), qfun, p = q90_lo)
from_FP_hi90B <- apply(boot_C_rw[, idx_FP, , ], c(1, 2), qfun, p = q90_hi)
from_FP_lo68B <- apply(boot_C_rw[, idx_FP, , ], c(1, 2), qfun, p = q68_lo)
from_FP_hi68B <- apply(boot_C_rw[, idx_FP, , ], c(1, 2), qfun, p = q68_hi)

# Truncar bandas al dominio [0, 100]
from_FP_lo90B <- pmax(0, from_FP_lo90B)
from_FP_hi90B <- pmin(100, from_FP_hi90B)
from_FP_lo68B <- pmax(0, from_FP_lo68B)
from_FP_hi68B <- pmin(100, from_FP_hi68B)

# orden_FP_con_own <- c("FP",         # own primero
#                       "BancoCentral", "Bancos", "Seguros", 
#                       "FMM", "OFIs", "FNM", "RestoMundo")

df_from_FP <- tibble(
  fecha    = rep(as.Date(fechas_rw), K),
  sector   = rep(nombre_vars, each = n_dates),
  BC       = as.vector(from_FP_obs),   # ahora es el observado, no el BC
  lo90     = as.vector(from_FP_lo90B), hi90 = as.vector(from_FP_hi90B),
  lo68     = as.vector(from_FP_lo68B), hi68 = as.vector(from_FP_hi68B)
) %>% 
  #filter(sector != "FP") %>%
  mutate(sector_label = ifelse(sector == "FP", "FP(own)", as.character(sector)),
         sector_label = factor(sector_label, levels = c("FP(own)", "BancoCentral", 
                                                        "Bancos", "Seguros",
                                                        "FMM", "OFIs", 
                                                        "FNM", "RestoMundo"))
         )

write_csv(df_from_FP, file.path(ruta_output, "bcch_21_pairwise_FROM_FP.csv"))

# --- (b) TO desde BancoCentral: C[i, BCCh] para i != BCCh, por fecha ----------
to_BCCh_obs  <- t(sapply(1:n_dates, function(t) CT_rw_obs[, idx_BCCh, t]))
colnames(to_BCCh_obs) <- nombre_vars

to_BCCh_bar  <- apply(boot_C_rw[, , idx_BCCh, ], c(1, 2), mean,   na.rm = TRUE)
to_BCCh_lo90 <- apply(boot_C_rw[, , idx_BCCh, ], c(1, 2), qfun, p = q90_hi)
to_BCCh_hi90 <- apply(boot_C_rw[, , idx_BCCh, ], c(1, 2), qfun, p = q90_lo)
to_BCCh_lo68 <- apply(boot_C_rw[, , idx_BCCh, ], c(1, 2), qfun, p = q68_hi)
to_BCCh_hi68 <- apply(boot_C_rw[, , idx_BCCh, ], c(1, 2), qfun, p = q68_lo)

to_BCCh_BC    <- 2 * to_BCCh_obs - to_BCCh_bar
to_BCCh_lo90B <- 2 * to_BCCh_obs - to_BCCh_lo90
to_BCCh_hi90B <- 2 * to_BCCh_obs - to_BCCh_hi90
to_BCCh_lo68B <- 2 * to_BCCh_obs - to_BCCh_lo68
to_BCCh_hi68B <- 2 * to_BCCh_obs - to_BCCh_hi68

to_BCCh_lo90B <- pmax(0, to_BCCh_lo90B)
to_BCCh_hi90B <- pmin(100, to_BCCh_hi90B)
to_BCCh_lo68B <- pmax(0, to_BCCh_lo68B)
to_BCCh_hi68B <- pmin(100, to_BCCh_hi68B)

df_to_BCCh <- tibble(
  fecha    = rep(as.Date(fechas_rw), K),
  sector   = rep(nombre_vars, each = n_dates),
  #obs      = as.vector(to_BCCh_obs),
  BC      = as.vector(to_BCCh_obs),
  #BC       = as.vector(to_BCCh_BC),
  lo90     = as.vector(to_BCCh_lo90B), hi90 = as.vector(to_BCCh_hi90B),
  lo68     = as.vector(to_BCCh_lo68B), hi68 = as.vector(to_BCCh_hi68B)
) %>% 
  filter(sector != "BancoCentral") %>%
  mutate(sector = factor(sector, levels = setdiff(nombre_vars, "BancoCentral")))

write_csv(df_to_BCCh, file.path(ruta_output, "bcch_22_pairwise_TO_BancoCentral.csv"))

# --- Graficos ----------------------------------------------------------------
# Figura 1: FROM hacia FP
orden_narrativo_FROM_FP <- c("BancoCentral", "Bancos", "FP(own)", "Seguros", 
                             "FMM", "FNM", "OFIs", "RestoMundo")
df_from_FP <- df_from_FP %>% 
  mutate(sector_label = factor(sector_label, levels = orden_narrativo_FROM_FP))

p_from_FP <- ggplot(df_from_FP, aes(x = fecha)) +
  #geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#4682B4", alpha = 0.18) +
  #geom_ribbon(aes(ymin = lo68, ymax = hi68), fill = "#4682B4", alpha = 0.35) +
  geom_line(aes(y = BC), color = "steelblue", linewidth = 0.7) +
  facet_wrap(~ sector_label, scales = "free_y", ncol = 2) +
  eventos() +
  scale_x_date(date_breaks = "3 years",
               labels = function(x) ifelse(is.na(x), "", lubridate::year(x))) +
  labs(title    = "Contribuciones recibidas por Fondos de Pensiones, por sector emisor",
       subtitle = "Pairwise FROM hacia FP, Chile 2011T3-2025T3",
       x = NULL, y = "% de varianza de FP explicada por shocks del sector emisor") +
  tema_g
ggsave(file.path(ruta_output, "figs", "bcch_20_pairwise_FROM_FP.png"),
       p_from_FP, width = 6.5, height = 7.5, dpi = 300, bg = "white")

# Figura 2: TO desde BancoCentral
orden_narrativo_TO_BC <- c("FP","Bancos", "Seguros", 
                     "FMM", "FNM", "OFIs", "RestoMundo")
df_to_BCCh <- df_to_BCCh %>% 
  mutate(sector = factor(sector, levels = orden_narrativo_TO_BC))

p_to_BCCh <- ggplot(df_to_BCCh, aes(x = fecha)) +
  #geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#4682B4", alpha = 0.18) +
  #geom_ribbon(aes(ymin = lo68, ymax = hi68), fill = "#4682B4", alpha = 0.35) +
  geom_line(aes(y = BC), color = "steelblue", linewidth = 0.7) +
  facet_wrap(~ sector, scales = "free_y", ncol = 2) +
  eventos() +
  scale_x_date(date_breaks = "3 years",
               labels = function(x) ifelse(is.na(x), "", lubridate::year(x))) +
  labs(title    = "Contribuciones del Banco Central a la varianza de otros sectores",
       subtitle = "Pairwise TO desde BancoCentral, Chile 2011T3-2025T3",
       x = NULL, y = "% de varianza del sector receptor explicada por shocks del Banco Central") +
  tema_g
ggsave(file.path(ruta_output, "figs", "bcch_21_pairwise_TO_BancoCentral.png"),
       p_to_BCCh, width = 6.5, height = 7.5, dpi = 300, bg = "white")

cat("  - bcch_20_pairwise_FROM_FP.png             FROM hacia FP, 7 facets\n")
cat("  - bcch_21_pairwise_TO_BancoCentral.png     TO desde BCCh, 7 facets\n")

# =============================================================================
# (6) FOREST PLOTS DE CONTRASTES PAIRWISE - ENFOQUE RETIROS
# =============================================================================
cat("\nGenerando contrastes pairwise para Retiros...\n")

# Definir contrastes filtrados a Retiros como referencia
contrastes_retiros <- tibble::tribble(
  ~A,              ~B,
  "Retiros",       "Pre-Estallido",
  "Retiros",       "Estallido",
  "Post-Retiros",  "Retiros"
)

# Strings de niveles para ordenar el factor correctamente
niveles_contrastes <- paste0(contrastes_retiros$A, " vs ", contrastes_retiros$B)

hacer_contraste <- function(boot_A, boot_B, obs_A, obs_B) {
  ic <- ic_basico_escalar(obs_A - obs_B, boot_A - boot_B)
  diff_dist <- boot_A - boot_B
  p_val <- 2 * min(mean(diff_dist <= 0, na.rm = TRUE),
                   mean(diff_dist >= 0, na.rm = TRUE))
  list(obs = obs_A - obs_B, BC = ic$theta_BC,
       ic68_lo = ic$lo68, ic68_hi = ic$hi68,
       ic90_lo = ic$lo90, ic90_hi = ic$hi90,
       p_value = p_val)
}

# Contrastes TCI - iterando sobre contrastes_retiros, NO contrastes
contrast_tci <- bind_rows(lapply(1:nrow(contrastes_retiros), function(k) {
  A <- contrastes_retiros$A[k]; B_ <- contrastes_retiros$B[k]
  ct <- hacer_contraste(tci_boot_reg[, A], tci_boot_reg[, B_],
                        tci_obs_reg[A], tci_obs_reg[B_])
  tibble(contraste = paste0(A, " vs ", B_), 
         metric = "TCI", sector = NA_character_,
         obs = ct$obs, BC = ct$BC,
         ic68_lo = ct$ic68_lo, ic68_hi = ct$ic68_hi,
         ic90_lo = ct$ic90_lo, ic90_hi = ct$ic90_hi,
         p_value = ct$p_value,
         sig_68 = !(ct$ic68_lo <= 0 & ct$ic68_hi >= 0),
         sig_90 = !(ct$ic90_lo <= 0 & ct$ic90_hi >= 0))
})) %>%
  mutate(contraste = factor(contraste, levels = niveles_contrastes))

# Contrastes NET por sector - iterando sobre contrastes_retiros
contrast_net <- bind_rows(lapply(1:nrow(contrastes_retiros), function(k) {
  A <- contrastes_retiros$A[k]; B_ <- contrastes_retiros$B[k]
  bind_rows(lapply(1:K, function(s) {
    sname <- nombre_vars[s]
    ct <- hacer_contraste(net_boot_reg[, s, A], net_boot_reg[, s, B_],
                          net_obs_reg[s, A], net_obs_reg[s, B_])
    tibble(contraste = paste0(A, " vs ", B_), 
           metric = "NET", sector = sname,
           obs = ct$obs, BC = ct$BC,
           ic68_lo = ct$ic68_lo, ic68_hi = ct$ic68_hi,
           ic90_lo = ct$ic90_lo, ic90_hi = ct$ic90_hi,
           p_value = ct$p_value,
           sig_68 = !(ct$ic68_lo <= 0 & ct$ic68_hi >= 0),
           sig_90 = !(ct$ic90_lo <= 0 & ct$ic90_hi >= 0))
  }))
})) %>%
  mutate(contraste = factor(contraste, levels = niveles_contrastes))

# Guardar CSV
df_contrastes <- bind_rows(contrast_tci, contrast_net) %>%
  mutate(stars = case_when(
    p_value < 0.01 ~ "***", p_value < 0.05 ~ "**",
    p_value < 0.10 ~ "*",   p_value < 0.32 ~ ".",   TRUE ~ ""
  ))
write_csv(df_contrastes, file.path(ruta_output, "bcch_20_regime_contrasts.csv"))

cat("\nContrastes TCI (con IC basico):\n")
print(contrast_tci %>%
        mutate(across(c(obs, BC, ic68_lo, ic68_hi, ic90_lo, ic90_hi), ~round(., 2)),
               p_value = round(p_value, 4)) %>%
        select(contraste, obs, BC, ic90_lo, ic90_hi, p_value, sig_90))

# --- Plot medias TCI por regimen ---
p_reg_tci <- ggplot(df_regime_tci, aes(x = regimen, y = BC)) +
  geom_col(fill = "steelblue", alpha = 0.4, width = 0.6) +
  geom_errorbar(aes(ymin = ic90_lo, ymax = ic90_hi), 
                width = 0.15, linewidth = 0.4, color = "steelblue4") +
  geom_errorbar(aes(ymin = ic68_lo, ymax = ic68_hi), 
                width = 0, linewidth = 2, color = "steelblue4") +
  geom_point(size = 2.5, color = "steelblue4") +
  geom_text(aes(label = sprintf("%.1f", BC)), 
            hjust = -0.3, vjust = -0.6, size = 2.8, 
            color = "gray20", family = "Times New Roman") +
  labs(title = "TCI medio por régimen",
       subtitle = sprintf("Bandas: 68%% (gruesa) y 90%% (fina). B=%d, W=%d", B_rolling, WINDOW), 
       x = NULL, y = "TCI (%)") +
  tema_g
ggsave(file.path(ruta_output, "figs", "bcch_17_regime_means_tci_BC.png"),
       p_reg_tci, width = 9, height = 5.5, dpi = 150, bg = "white")

# --- Forest plot TCI ---
p_forest_tci <- ggplot(contrast_tci,
                       aes(x = BC, y = fct_rev(contraste))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = ic90_lo, xmax = ic90_hi, color = sig_90), 
                 height = 0, linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ic68_lo, xmax = ic68_hi, color = sig_90), 
                 height = 0, linewidth = 2) +
  geom_point(aes(color = sig_90), size = 2.8) +
  geom_text(aes(label = sprintf("%.1f", BC)), 
            hjust = -0.3, vjust = -0.6, size = 2.8, 
            color = "gray20", family = "Times New Roman") +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "steelblue4")) +
  labs(title = "Diferencias en TCI entre regímenes",
       subtitle = "Diferencia A − B en puntos porcentuales. Bandas: 68% (gruesa) y 90% (fina).",
       x = "Diferencia (A − B), pp", y = NULL) +
  tema_g +
  theme(legend.position = "none")

ggsave(file.path(ruta_output, "figs", "bcch_22_contrastes_TCI.png"),
       p_forest_tci, width = 9, height = 5, dpi = 150, bg = "white")

# --- Forest plot NET por sector ---
p_forest_net <- ggplot(contrast_net,
                       aes(x = BC, y = fct_rev(factor(sector, levels = nombre_vars)))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = ic90_lo, xmax = ic90_hi, color = sig_90), 
                 height = 0, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ic68_lo, xmax = ic68_hi, color = sig_90), 
                 height = 0, linewidth = 1.5) +
  geom_point(aes(color = sig_90), size = 2) +
  facet_wrap(~ contraste, ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "steelblue4")) +
  labs(title = "Diferencias en connectividad neta entre regímenes, por sector",
       subtitle = "Diferencia A − B en puntos porcentuales. Bandas: 68% (gruesa) y 90% (fina).",
       x = "Diferencia (A − B), pp", y = NULL) +
  tema_g +
  theme(legend.position = "none")

ggsave(file.path(ruta_output, "figs", "bcch_23_contrastes_NET.png"),
       p_forest_net, width = 9, height = 9, dpi = 150, bg = "white")

cat("  - bcch_22_contrastes_TCI.png     Forest plot diferencias TCI\n")
cat("  - bcch_23_contrastes_NET.png     Forest plot diferencias NET por sector\n")


# =============================================================================
# NOTA SOBRE EXTENSIONES OPCIONALES (no implementadas, mencionadas para tesis):
# =============================================================================
# 1. Correccion de sesgo de Kilian (1998): re-estimar el VAR con coeficientes
#    bias-corrected antes de bootstrap. Util cuando T es chico (aqui T~89).
# 2. Bootstrap-after-bootstrap (Greenwood-Nimmo et al. 2024 Sec. 3.2): para
#    coverage nominal correcta del IC. Costo computacional: B^2 reps.
# 3. Block bootstrap (Politis-Romano) si los residuos siguen mostrando algun
#    grado de dependencia serial mas alla del VAR(P_LAG).
# 4. Paralelizacion: envolver el loop con future.apply::future_lapply()
#    despues de plan(multisession, workers = parallel::detectCores() - 1).
# =============================================================================