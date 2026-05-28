# =============================================================================
# Script: 02_bootstrap_dy_BCCh.R   (v5: bandas duales + regimenes + bias corr.)
# Autor:  Carlos Gonzalez | Mayo 2026
#
# CAMBIO PRINCIPAL v5 vs v4:
# Aplica IC bootstrap basico (Hall 1992) en vez de percentil, lo que corrige
# el sesgo de muestra pequena visible en las salidas previas (la banda
# percentil quedaba sistematicamente por encima del estimado puntual del TCI
# en el periodo Pre-Estallido). El sesgo viene del bias de OLS en VAR (Nickell
# 1981) propagado a traves del FEVD: A_hat subestima la persistencia, lo que
# infla artificialmente el TCI tanto en el estimado como en el bootstrap.
#
# Formulas usadas:
#   bias_hat   = mean(theta*) - theta_hat                  (sesgo estimado)
#   theta_BC   = theta_hat - bias_hat = 2*theta_hat - mean(theta*)
#   IC_basico  = [2*theta_hat - q_hi(theta*), 2*theta_hat - q_lo(theta*)]
#
# El IC basico es matematicamente equivalente al percentil clasico bajo el
# supuesto pivotal estandar, pero corrige automaticamente el sesgo. Es lo
# recomendado en Davison & Hinkley (1997) Cap. 5 para estadisticos sesgados.
#
# Una alternativa mas rigurosa seria la correccion de Kilian (1998) sobre A
# (bootstrap-after-bootstrap), pero requiere ~50% mas de computo y para una
# tesis de magister el basico es suficiente y se interpreta igual.
#
# REFERENCIAS:
#   Davison & Hinkley (1997). Bootstrap Methods and Their Application. Cap. 5.
#   Hall (1992). The Bootstrap and Edgeworth Expansion. Sec. 3.5.
#   Kilian (1998). Small-sample confidence intervals for impulse response
#     functions. Review of Economics and Statistics 80, 218-230.
#   Greenwood-Nimmo, Kocenda, Nguyen (2024). Detecting statistically significant
#     changes in connectedness. Economic Modelling 140, 106843.
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
ruta_output_1 <- ruta/"3.2_resultado_dy/salidas_1"

# =============================================================================
# PASO 1 - SETUP DEL ENTORNO
# =============================================================================
paquetes <- c("ConnectednessApproach", "vars", "zoo",
              "tidyverse", "lubridate", "progress")
paquetes_faltantes <- paquetes[!paquetes %in% installed.packages()[, "Package"]]
if (length(paquetes_faltantes) > 0) install.packages(paquetes_faltantes)

library(ConnectednessApproach)
library(vars)
library(zoo)
library(tidyverse)
library(lubridate)
library(progress)

dir.create(ruta_output,                           showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ruta_output, "figs"),        showWarnings = FALSE)

options(scipen = 999)
set.seed(2024)

# =============================================================================
# PARAMETROS
# =============================================================================
B_static   <- 2000
B_rolling  <- 1000
NFORE      <- 8
WINDOW     <- 34
ALPHA_90   <- 0.10
ALPHA_68   <- 0.32
BOOT_TYPE  <- "wild"

# Regimenes (post-2011 porque rolling-window arranca en 2011T3)
regimenes <- tibble::tribble(
  ~regimen,         ~inicio,         ~fin,
  "Pre-Estallido",  "2011-09-01",    "2019-08-30",   # ~33 obs
  "Estallido",      "2019-09-01",    "2020-03-30",   # ~3  obs
  "Retiros",        "2020-03-01",    "2021-09-30",   # ~6  obs (considero trimestre inmediatamente anterior y posterior)
  "Post-Retiros",   "2022-01-01",    "2025-09-30"    # ~15 obs
) %>% mutate(inicio = as.Date(inicio), fin = as.Date(fin))

contrastes <- tibble::tribble(
  ~A,              ~B,
  "Estallido",     "Pre-Estallido",
  "Retiros",       "Pre-Estallido",
  "Retiros",       "Estallido",
  "Post-Retiros",  "Pre-Estallido",
  "Post-Retiros",  "Retiros"
)

cat("\n===============================================================\n")
cat("BOOTSTRAP DY v5 - IC basico (corr. sesgo) + bandas duales + regimenes\n")
cat("===============================================================\n")

# =============================================================================
# DATA
# =============================================================================
dy_wide <- read_csv(file.path(ruta_output_1, "01_serie_var_wide.csv"),
                    show_col_types = FALSE)

orden_sectores <- c("S129", "S122", "S2", "S128", "S123", "S124", "S125_S126", "S121")
nombres_sectores <- c(
  "S129"="FP", "S122"="Bancos", "S2"="RestoMundo", "S128"="Seguros",
  "S123"="FMM", "S124"="FNM", "S125_S126"="OFIs", "S121"="BancoCentral"
)

periodo_a_fecha <- function(periodo) {
  year <- as.numeric(substr(periodo, 1, 4))
  trim <- as.numeric(substr(periodo, 6, 6))
  mes_inicio <- (trim - 1) * 3 + 1
  as.Date(paste0(year, "-", sprintf("%02d", mes_inicio), "-01"))
}

fechas   <- periodo_a_fecha(dy_wide$periodo)
Y_matrix <- as.matrix(dy_wide[, orden_sectores])
colnames(Y_matrix) <- nombres_sectores[orden_sectores]
Y_zoo    <- zoo(Y_matrix, order.by = fechas)

P_LAG <- as.integer(VARselect(coredata(Y_zoo), lag.max = 4,
                              type = "const")$selection["SC(n)"])

q90_lo <- ALPHA_90 / 2; q90_hi <- 1 - ALPHA_90 / 2
q68_lo <- ALPHA_68 / 2; q68_hi <- 1 - ALPHA_68 / 2

cat(sprintf("VAR p=%d | nfore=%d | W=%d | B_static=%d | B_rolling=%d\n",
            P_LAG, NFORE, WINDOW, B_static, B_rolling))

# =============================================================================
# HELPERS
# =============================================================================
extraer_CT_estatica <- function(out_obj, K) {
  C <- NULL
  for (slot in c("CT", "gFEVD", "TOTAL", "FEVD", "PHI")) {
    if (!is.null(out_obj[[slot]])) {
      cand <- out_obj[[slot]]
      if (is.array(cand) && length(dim(cand)) == 3) {
        if (dim(cand)[3] == 1) cand <- cand[, , 1] else next
      }
      cand <- tryCatch(as.matrix(cand), error = function(e) NULL)
      if (!is.null(cand) && is.numeric(cand) && nrow(cand) == K && ncol(cand) == K) {
        C <- cand; break
      }
    }
  }
  if (is.null(C)) stop("No se extrajo matriz K x K")
  storage.mode(C) <- "double"
  if (max(abs(C), na.rm = TRUE) <= 1.5) C <- C * 100
  C
}

extraer_CT_rolling <- function(out_obj, K) {
  CT <- NULL
  for (slot in c("CT", "gFEVD", "TOTAL", "FEVD")) {
    if (!is.null(out_obj[[slot]])) {
      cand <- out_obj[[slot]]
      if (is.array(cand) && length(dim(cand)) == 3 &&
          dim(cand)[1] == K && dim(cand)[2] == K) {
        CT <- cand; break
      }
    }
  }
  if (is.null(CT)) stop("No se extrajo CT rolling")
  storage.mode(CT) <- "double"
  if (max(abs(CT), na.rm = TRUE) <= 1.5) CT <- CT * 100
  CT
}

calc_agregados <- function(C) {
  d <- diag(C)
  from_i <- rowSums(C) - d
  to_j   <- colSums(C) - d
  list(from = from_i, to = to_j, net = to_j - from_i, tci = mean(from_i))
}

simular_Y_estrella <- function(Y, p, Ahat_list, nu_hat, resids, type = "wild") {
  Tt <- nrow(Y); k <- ncol(Y)
  Yst <- matrix(NA_real_, Tt, k); colnames(Yst) <- colnames(Y)
  Yst[1:p, ] <- as.matrix(Y[1:p, , drop = FALSE])
  Tres <- nrow(resids)
  if (type == "wild") {
    a <- -(sqrt(5)-1)/2; bb <- (sqrt(5)+1)/2
    pa <- (sqrt(5)+1)/(2*sqrt(5))
    w <- sample(c(a, bb), Tres, replace = TRUE, prob = c(pa, 1-pa))
    eps <- resids * w
  } else {
    eps <- resids[sample.int(Tres, Tres, replace = TRUE), , drop = FALSE]
  }
  for (t in (p+1):Tt) {
    yt <- nu_hat
    for (lag in 1:p) yt <- yt + Ahat_list[[lag]] %*% Yst[t-lag, ]
    yt <- yt + eps[t-p, ]
    Yst[t, ] <- as.numeric(yt)
  }
  Yst
}

qfun <- function(x, p) quantile(x, probs = p, na.rm = TRUE, names = FALSE)

# Funcion CENTRAL: dado theta_hat (escalar o vector) y matriz boot (replicas en
# columnas o filas segun caso), devuelve la lista con:
#   theta_bar : media bootstrap
#   bias      : theta_bar - theta_hat (sesgo estimado)
#   theta_BC  : theta_hat - bias = 2*theta_hat - theta_bar
#   lo90/hi90, lo68/hi68 : IC basico (reflejado, corrige sesgo)
# Cuando theta_hat es un vector y boot es matriz [n x B], opera por fila.
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

# =============================================================================
# AJUSTE VAR
# =============================================================================
fit_full    <- VAR(coredata(Y_zoo), p = P_LAG, type = "const")
Ahat_full   <- Acoef(fit_full)
nu_full     <- sapply(fit_full$varresult, function(eq) unname(coef(eq)["const"]))
res_full    <- residuals(fit_full)
K           <- ncol(Y_zoo)
nombre_vars <- colnames(Y_zoo)

# =============================================================================
# (1) BOOTSTRAP ESTATICO
# =============================================================================
cat("\n--- (1) BOOTSTRAP ESTATICO ---\n")

dy_full <- ConnectednessApproach(
  x = Y_zoo, nlag = P_LAG, nfore = NFORE,
  model = "VAR", connectedness = "Time", window.size = NULL
)
C_obs <- extraer_CT_estatica(dy_full, K)
rownames(C_obs) <- colnames(C_obs) <- nombre_vars
agg_obs <- calc_agregados(C_obs)

boot_fevd <- array(NA_real_, c(K, K, B_static), dimnames = list(nombre_vars, nombre_vars, NULL))
boot_to   <- matrix(NA_real_, K, B_static, dimnames = list(nombre_vars, NULL))
boot_from <- matrix(NA_real_, K, B_static, dimnames = list(nombre_vars, NULL))
boot_net  <- matrix(NA_real_, K, B_static, dimnames = list(nombre_vars, NULL))
boot_tci  <- rep(NA_real_, B_static)

pb <- progress_bar$new(total = B_static, format = "[:bar] :percent eta :eta", width = 60)
for (b in seq_len(B_static)) {
  pb$tick()
  Yst <- simular_Y_estrella(coredata(Y_zoo), P_LAG, Ahat_full, nu_full, res_full, BOOT_TYPE)
  out <- tryCatch(ConnectednessApproach(
    x = zoo(Yst, order.by = index(Y_zoo)), nlag = P_LAG, nfore = NFORE,
    model = "VAR", connectedness = "Time", window.size = NULL
  ), error = function(e) NULL)
  if (is.null(out)) next
  Cb <- tryCatch(extraer_CT_estatica(out, K), error = function(e) NULL)
  if (is.null(Cb)) next
  ag <- calc_agregados(Cb)
  boot_fevd[, , b] <- Cb
  boot_to[, b]   <- ag$to
  boot_from[, b] <- ag$from
  boot_net[, b]  <- ag$net
  boot_tci[b]    <- ag$tci
}

# IC basico para cada celda y cada agregado
# Para FEVD K x K:
fevd_bar  <- apply(boot_fevd, c(1,2), mean, na.rm = TRUE)
fevd_lo90 <- 2*C_obs - apply(boot_fevd, c(1,2), qfun, p = q90_hi)
fevd_hi90 <- 2*C_obs - apply(boot_fevd, c(1,2), qfun, p = q90_lo)
fevd_BC   <- 2*C_obs - fevd_bar
dimnames(fevd_BC) <- dimnames(fevd_lo90) <- dimnames(fevd_hi90) <- dimnames(C_obs)

# Para vectores TO, FROM, NET y escalar TCI
ic_to   <- lapply(1:K, function(s) ic_basico_escalar(agg_obs$to[s],   boot_to[s, ]))
ic_from <- lapply(1:K, function(s) ic_basico_escalar(agg_obs$from[s], boot_from[s, ]))
ic_net  <- lapply(1:K, function(s) ic_basico_escalar(agg_obs$net[s],  boot_net[s, ]))
ic_tci  <- ic_basico_escalar(agg_obs$tci, boot_tci)

# Diagnostico de sesgo
cat(sprintf("TCI puntual (theta_hat) : %.2f%%\n", agg_obs$tci))
cat(sprintf("TCI media bootstrap     : %.2f%%\n", ic_tci$theta_bar))
cat(sprintf("Sesgo estimado          : %+.2f pp\n", ic_tci$bias))
cat(sprintf("TCI bias-corrected      : %.2f%%\n", ic_tci$theta_BC))
cat(sprintf("IC 90%% basico (corr.)   : [%.2f, %.2f]\n", ic_tci$lo90, ic_tci$hi90))

# Tabla "X.XX [lo, hi]" con bias-corrected en el centro
fmt_celda_BC <- function(BC, lo, hi, d = 2) sprintf("%.*f [%.*f, %.*f]", d, BC, d, lo, d, hi)

tabla_principal <- matrix("", K, K+1, dimnames = list(nombre_vars, c(nombre_vars, "FROM")))
for (i in 1:K) {
  for (j in 1:K) {
    tabla_principal[i,j] <- fmt_celda_BC(fevd_BC[i,j], fevd_lo90[i,j], fevd_hi90[i,j])
  }
  tabla_principal[i, K+1] <- fmt_celda_BC(ic_from[[i]]$theta_BC, ic_from[[i]]$lo90, ic_from[[i]]$hi90)
}
fila_to  <- c(sapply(1:K, function(j) fmt_celda_BC(ic_to[[j]]$theta_BC,  ic_to[[j]]$lo90,  ic_to[[j]]$hi90)),  "")
fila_net <- c(sapply(1:K, function(j) fmt_celda_BC(ic_net[[j]]$theta_BC, ic_net[[j]]$lo90, ic_net[[j]]$hi90)), "")
fila_tci <- c(rep("", K-1), "TCI:", fmt_celda_BC(ic_tci$theta_BC, ic_tci$lo90, ic_tci$hi90))
tabla_full <- rbind(tabla_principal, TO = fila_to, NET = fila_net, TCI = fila_tci)

write.csv(tabla_full, file.path(ruta_output, "bcch_15_bootstrap_tabla_DY_IC.csv"), row.names = TRUE)
write.csv(C_obs,     file.path(ruta_output, "bcch_15a_fevd_punto.csv"),       row.names = TRUE)
write.csv(fevd_BC,   file.path(ruta_output, "bcch_15a2_fevd_bias_corrected.csv"), row.names = TRUE)
write.csv(fevd_lo90, file.path(ruta_output, "bcch_15b_fevd_lower90.csv"),     row.names = TRUE)
write.csv(fevd_hi90, file.path(ruta_output, "bcch_15c_fevd_upper90.csv"),     row.names = TRUE)
write_csv(tibble(
  sector  = nombre_vars,
  TO      = agg_obs$to,   TO_BC   = sapply(ic_to,   `[[`, "theta_BC"),
  TO_lo90 = sapply(ic_to, `[[`, "lo90"),    TO_hi90 = sapply(ic_to, `[[`, "hi90"),
  FROM    = agg_obs$from, FROM_BC = sapply(ic_from, `[[`, "theta_BC"),
  FROM_lo90 = sapply(ic_from, `[[`, "lo90"), FROM_hi90 = sapply(ic_from, `[[`, "hi90"),
  NET     = agg_obs$net,  NET_BC  = sapply(ic_net,  `[[`, "theta_BC"),
  NET_lo90 = sapply(ic_net, `[[`, "lo90"), NET_hi90 = sapply(ic_net, `[[`, "hi90")
), file.path(ruta_output, "bcch_15d_directional_IC.csv"))

# =============================================================================
# (2-3) BOOTSTRAP ROLLING
# =============================================================================
cat("\n--- (2-3) BOOTSTRAP ROLLING ---\n")

dy_rw_obs <- ConnectednessApproach(
  x = Y_zoo, nlag = P_LAG, nfore = NFORE,
  model = "VAR", connectedness = "Time", window.size = WINDOW
)
CT_rw_obs <- extraer_CT_rolling(dy_rw_obs, K)
n_dates   <- dim(CT_rw_obs)[3]
fechas_rw <- tail(index(Y_zoo), n_dates)

TCI_obs_rw <- numeric(n_dates)
NET_obs_rw <- matrix(NA_real_, n_dates, K, dimnames = list(NULL, nombre_vars))
for (t in 1:n_dates) {
  ag <- calc_agregados(CT_rw_obs[, , t])
  TCI_obs_rw[t]   <- ag$tci
  NET_obs_rw[t, ] <- ag$net
}
cat(sprintf("Rolling: %d fechas, %s a %s\n", n_dates, format(min(fechas_rw)), format(max(fechas_rw))))

boot_tci_rw <- matrix(NA_real_, n_dates, B_rolling)
boot_net_rw <- array(NA_real_, c(n_dates, K, B_rolling),
                     dimnames = list(NULL, nombre_vars, NULL))

# Nuevo: array completo de matrices C por replica (para pairwise: FP FROM y BC TO)
boot_C_rw <- array(NA_real_, c(n_dates, K, K, B_rolling),
                   dimnames = list(NULL, nombre_vars, nombre_vars, NULL))

pb2 <- progress_bar$new(total = B_rolling, format = "[:bar] :percent eta :eta", width = 60)
for (b in seq_len(B_rolling)) {
  pb2$tick()
  Yst <- simular_Y_estrella(coredata(Y_zoo), P_LAG, Ahat_full, nu_full, res_full, BOOT_TYPE)
  out <- tryCatch(ConnectednessApproach(
    x = zoo(Yst, order.by = index(Y_zoo)), nlag = P_LAG, nfore = NFORE,
    model = "VAR", connectedness = "Time", window.size = WINDOW
  ), error = function(e) NULL)
  if (is.null(out)) next
  CT_b <- tryCatch(extraer_CT_rolling(out, K), error = function(e) NULL)
  if (is.null(CT_b) || dim(CT_b)[3] != n_dates) next
  for (t in 1:n_dates) {
    ag <- calc_agregados(CT_b[, , t])
    boot_tci_rw[t, b]   <- ag$tci
    boot_net_rw[t, , b] <- ag$net
    boot_C_rw[t, , , b]  <- CT_b[, , t]   # nuevo: guardar matriz pairwise
    
  }
}

# IC basico para TCI rolling (por fecha)
tci_bar_rw  <- rowMeans(boot_tci_rw, na.rm = TRUE)
tci_bias_rw <- tci_bar_rw - TCI_obs_rw
tci_BC_rw   <- TCI_obs_rw - tci_bias_rw       # = 2*TCI_obs_rw - tci_bar_rw
tci_lo90_BC <- 2*TCI_obs_rw - apply(boot_tci_rw, 1, qfun, p = q90_hi)
tci_hi90_BC <- 2*TCI_obs_rw - apply(boot_tci_rw, 1, qfun, p = q90_lo)
tci_lo68_BC <- 2*TCI_obs_rw - apply(boot_tci_rw, 1, qfun, p = q68_hi)
tci_hi68_BC <- 2*TCI_obs_rw - apply(boot_tci_rw, 1, qfun, p = q68_lo)

cat(sprintf("Sesgo medio TCI rolling: %+.2f pp (W=34 -> sesgo mas grande que estatico)\n",
            mean(tci_bias_rw, na.rm = TRUE)))

df_tci_boot <- tibble(
  fecha    = as.Date(fechas_rw),
  tci_obs  = TCI_obs_rw,             # estimado puntual (uncorrected, mantener para comparacion)
  tci_bar  = tci_bar_rw,             # media bootstrap
  tci_BC   = tci_BC_rw,              # bias-corrected (esto va en el grafico)
  tci_lo90 = tci_lo90_BC, tci_hi90 = tci_hi90_BC,
  tci_lo68 = tci_lo68_BC, tci_hi68 = tci_hi68_BC
)
write_csv(df_tci_boot, file.path(ruta_output, "bcch_16_bootstrap_tci_rolling.csv"))

# IC basico para NET rolling (por fecha y sector)
net_bar_rw  <- apply(boot_net_rw, c(1,2), mean, na.rm = TRUE)
net_BC_rw   <- 2*NET_obs_rw - net_bar_rw
net_lo90_BC <- 2*NET_obs_rw - apply(boot_net_rw, c(1,2), qfun, p = q90_hi)
net_hi90_BC <- 2*NET_obs_rw - apply(boot_net_rw, c(1,2), qfun, p = q90_lo)
net_lo68_BC <- 2*NET_obs_rw - apply(boot_net_rw, c(1,2), qfun, p = q68_hi)
net_hi68_BC <- 2*NET_obs_rw - apply(boot_net_rw, c(1,2), qfun, p = q68_lo)

df_net_boot <- tibble(
  fecha   = rep(as.Date(fechas_rw), K),
  sector  = rep(nombre_vars, each = n_dates),
  obs     = as.vector(NET_obs_rw),
  bar     = as.vector(net_bar_rw),
  BC      = as.vector(net_BC_rw),
  lo90    = as.vector(net_lo90_BC), hi90 = as.vector(net_hi90_BC),
  lo68    = as.vector(net_lo68_BC), hi68 = as.vector(net_hi68_BC)
) %>% mutate(sector = factor(sector, levels = nombre_vars))
write_csv(df_net_boot, file.path(ruta_output, "bcch_17_bootstrap_net_rolling.csv"))


saveRDS(list(
  # Arrays bootstrap principales
  boot_tci_rw = boot_tci_rw,
  boot_net_rw = boot_net_rw,
  boot_C_rw   = boot_C_rw,           # NUEVO: pairwise completo
  df_net_boot = df_net_boot,
  df_tci_boot = df_tci_boot,
  
  # Observados
  TCI_obs_rw  = TCI_obs_rw,
  NET_obs_rw  = NET_obs_rw,
  CT_rw_obs   = CT_rw_obs,           # NUEVO: matrices C observadas
  
  # Metadata
  fechas_rw   = fechas_rw,
  nombre_vars = nombre_vars,
  B_rolling   = B_rolling,
  WINDOW      = WINDOW,
  P_LAG       = P_LAG,
  NFORE       = NFORE,
  ALPHA_90    = ALPHA_90,
  ALPHA_68    = ALPHA_68,
  K           = K
),
file.path(ruta_output, "bcch_bootstrap_rolling.rds"))

# =============================================================================
# GRAFICOS
# =============================================================================
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
  "Pre-Estallido",  "2011-09-01",    "2019-09-30",   # ~33 obs
  "Estallido",      "2019-10-01",    "2020-06-30",   # ~3  obs
  "Retiros",        "2020-07-01",    "2021-12-31",   # ~6  obs
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
from_FP_obs  <- t(sapply(1:n_dates, function(t) CT_rw_obs[idx_FP, , t]))
colnames(from_FP_obs) <- nombre_vars

# Bootstrap: boot_C_rw[t, FP, j, b]
# Para cada (t, j) computar media y cuantiles bootstrap sobre b
from_FP_bar  <- apply(boot_C_rw[, idx_FP, , ], c(1, 2), mean,   na.rm = TRUE)
from_FP_lo90 <- apply(boot_C_rw[, idx_FP, , ], c(1, 2), qfun, p = q90_hi)
from_FP_hi90 <- apply(boot_C_rw[, idx_FP, , ], c(1, 2), qfun, p = q90_lo)
from_FP_lo68 <- apply(boot_C_rw[, idx_FP, , ], c(1, 2), qfun, p = q68_hi)
from_FP_hi68 <- apply(boot_C_rw[, idx_FP, , ], c(1, 2), qfun, p = q68_lo)

# IC basico bias-corrected
from_FP_BC    <- 2 * from_FP_obs - from_FP_bar
from_FP_lo90B <- 2 * from_FP_obs - from_FP_lo90
from_FP_hi90B <- 2 * from_FP_obs - from_FP_hi90
from_FP_lo68B <- 2 * from_FP_obs - from_FP_lo68
from_FP_hi68B <- 2 * from_FP_obs - from_FP_hi68


# Armar tibble long-format, incluyendo el sector FP (diagonal)
orden_FP_con_own <- c("BancoCentral", "Bancos", "Seguros", 
                      "FMM", "OFIs", "FNM", "RestoMundo",
                      "FP") # own incluido ahora

df_from_FP <- tibble(
  fecha    = rep(as.Date(fechas_rw), K),
  sector   = rep(nombre_vars, each = n_dates),
  obs      = as.vector(from_FP_obs),
  BC       = as.vector(from_FP_BC),
  lo90     = as.vector(from_FP_lo90B), hi90 = as.vector(from_FP_hi90B),
  lo68     = as.vector(from_FP_lo68B), hi68 = as.vector(from_FP_hi68B)
) %>% 
  #filter(sector != "FP") %>%
  mutate(sector = factor(sector, levels = orden_FP_con_own), 
         sector_label = ifelse(sector == "FP", "FP (own)", as.character(sector)))

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

df_to_BCCh <- tibble(
  fecha    = rep(as.Date(fechas_rw), K),
  sector   = rep(nombre_vars, each = n_dates),
  obs      = as.vector(to_BCCh_obs),
  BC       = as.vector(to_BCCh_BC),
  lo90     = as.vector(to_BCCh_lo90B), hi90 = as.vector(to_BCCh_hi90B),
  lo68     = as.vector(to_BCCh_lo68B), hi68 = as.vector(to_BCCh_hi68B)
) %>% 
  filter(sector != "BancoCentral") %>%
  mutate(sector = factor(sector, levels = setdiff(nombre_vars, "BancoCentral")))

write_csv(df_to_BCCh, file.path(ruta_output, "bcch_22_pairwise_TO_BancoCentral.csv"))

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

# --- Graficos ----------------------------------------------------------------
# Figura 1: FROM hacia FP
p_from_FP <- ggplot(df_from_FP, aes(x = fecha)) +
  geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#4682B4", alpha = 0.18) +
  geom_ribbon(aes(ymin = lo68, ymax = hi68), fill = "#4682B4", alpha = 0.35) +
  geom_line(aes(y = BC), color = "#1F3A60", linewidth = 0.7) +
  facet_wrap(~ sector, scales = "free_y", ncol = 2) +
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
p_to_BCCh <- ggplot(df_to_BCCh, aes(x = fecha)) +
  geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#A0826D", alpha = 0.18) +
  geom_ribbon(aes(ymin = lo68, ymax = hi68), fill = "#A0826D", alpha = 0.35) +
  geom_line(aes(y = BC), color = "#5C4033", linewidth = 0.7) +
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