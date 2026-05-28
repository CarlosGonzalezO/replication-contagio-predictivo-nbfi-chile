# =============================================================================
# TESIS MAGISTER ECONOMIA UC
# Intermediacion financiera no bancaria e interconexiones sectoriales en Chile
# =============================================================================
# Script: 01_dy_adaptive_lasso_bcch.R
# Autor:  Carlos Gonzalez
# Fecha:  Mayo 2026
#
# OBJETIVO:
# Replicar el pipeline DY (Pasos 1-5) usando Adaptive LASSO-VAR como
# estimador principal, en reemplazo del VAR irrestricto OLS.
#
# ESTRATEGIA DE ESTIMACION:
# El paquete ConnectednessApproach no implementa adaptive LASSO nativo.
# Se implementa manualmente en dos etapas por ecuacion del VAR:
#
#   Etapa 1 (OLS): estimar VAR irrestricto para obtener coeficientes iniciales
#                  hat_phi_OLS como pesos adaptativos.
#   Etapa 2 (aLASSO): estimar cada ecuacion del VAR con glmnet usando
#                     penalty.factor = 1 / |hat_phi_OLS| (adaptive weights).
#                     Seleccion de lambda por cross-validation (cv.glmnet).
#
# Luego se reconstruye la matriz Phi_hat a partir de los coeficientes aLASSO
# y se calcula la FEVD generalizada manualmente siguiendo DY (2012),
# generando los mismos indices (TCI, NET, TO, FROM, pairwise) que el VAR OLS.
#
# OPCION ALTERNATIVA (mas simple, menos control):
# Usar model = "LASSO" en ConnectednessApproach, que implementa LASSO
# estandar via BigVAR (sin adaptive weights). Ver seccion al final del script.
#
# INPUTS:
#   panel_var_dy.csv
#
# OUTPUTS (en /salidas_lasso/):
#   Mismos CSVs y graficos que el pipeline VAR original, prefijados con "lasso_"
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
ruta_output <- ruta / "3.2_resultado_dy/salidas_alasso"


# =============================================================================
# PASO 1 - SETUP
# =============================================================================

paquetes <- c(
  "ConnectednessApproach", "vars", "urca", "zoo",
  "tidyverse", "lubridate", "glmnet", "Matrix", "writexl"
)
paquetes_faltantes <- paquetes[!paquetes %in% installed.packages()[, "Package"]]
if (length(paquetes_faltantes) > 0) install.packages(paquetes_faltantes)

library(ConnectednessApproach)
library(vars)
library(urca)
library(zoo)
library(tidyverse)
library(lubridate)
library(glmnet)    # implementacion del adaptive LASSO
library(Matrix)
library(writexl)

dir.create(ruta_output,                      showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ruta_output, "figs"),   showWarnings = FALSE)

options(scipen = 999)
set.seed(42)

cat("\n===============================================================\n")
cat("DIEBOLD-YILMAZ: ADAPTIVE LASSO-VAR\n")
cat("===============================================================\n")


# =============================================================================
# PASO 2 - CARGA Y VALIDACION DE DATA
# (identico al script VAR original)
# =============================================================================
cat("\n--- PASO 2: CARGA Y VALIDACION ---\n")

archivo <- file.path(ruta_input, "panel_var_dy.csv")
if (!file.exists(archivo)) stop("Archivo no encontrado. Verifica pipeline Python.")

dy_long <- read_csv(archivo, show_col_types = FALSE)

dy_wide <- dy_long %>%
  select(periodo, sector_codigo, d_log_A_real) %>%
  pivot_wider(names_from = sector_codigo, values_from = d_log_A_real) %>%
  arrange(periodo) %>%
  filter(periodo != "2003T1")

write_csv(dy_wide, file.path(ruta_output, "lasso_bcch_01_serie_var_wide.csv"))

orden_sectores <- c("S129","S122", "S2", "S128","S123", "S124", "S125_S126","S121")
nombres_sectores <- c(
  "S129"      = "FP",
  "S122"      = "Bancos",
  "S2"        = "RestoMundo",
  "S128"      = "Seguros",
  "S123"      = "FMM",
  "S124"      = "FNM",
  "S125_S126" = "OFIs",
  "S121"      = "BancoCentral"
)

periodo_a_fecha <- function(periodo) {
  year  <- as.numeric(substr(periodo, 1, 4))
  trim  <- as.numeric(substr(periodo, 6, 6))
  mes   <- (trim - 1) * 3 + 1
  as.Date(paste0(year, "-", sprintf("%02d", mes), "-01"))
}

fechas   <- periodo_a_fecha(dy_wide$periodo)
Y_matrix <- as.matrix(dy_wide[, orden_sectores])
colnames(Y_matrix) <- nombres_sectores[orden_sectores]
Y_zoo    <- zoo(Y_matrix, order.by = fechas)


# =============================================================================
# PASO 3 - TESTS DE ESTACIONARIEDAD
# (identico al script VAR original)
# =============================================================================
cat("\n--- PASO 3: TESTS DE ESTACIONARIEDAD ---\n")

testear_estacionariedad <- function(serie, nombre_serie) {
  serie_limpia     <- na.omit(serie)
  adf              <- ur.df(serie_limpia, type = "drift", selectlags = "AIC")
  adf_rechaza_H0   <- adf@teststat[1] < adf@cval[1, "5pct"]
  kpss             <- ur.kpss(serie_limpia, type = "mu", lags = "short")
  kpss_rechaza_H0  <- kpss@teststat > kpss@cval[1, "5pct"]
  conclusion <- if (adf_rechaza_H0 & !kpss_rechaza_H0) "ESTACIONARIA"
  else if (!adf_rechaza_H0 & kpss_rechaza_H0) "NO ESTACIONARIA"
  else "AMBIGUA"
  data.frame(
    sector_codigo = nombre_serie,
    sector_nombre = nombres_sectores[nombre_serie],
    conclusion    = conclusion,
    stringsAsFactors = FALSE
  )
}

resultados_estacionariedad <- map_dfr(orden_sectores, function(s) {
  testear_estacionariedad(coredata(Y_zoo)[, nombres_sectores[s]], s)
})
write_csv(resultados_estacionariedad,
          file.path(ruta_output, "lasso_bcch_02_tests_estacionariedad.csv"))
print(resultados_estacionariedad)


# =============================================================================
# NUCLEO: FUNCIONES ADAPTIVE LASSO-VAR + FEVD GENERALIZADA
# =============================================================================

# -----------------------------------------------------------------------------
# estimar_alasso_var():
# Estima un VAR(p) ecuacion por ecuacion con adaptive LASSO.
# Devuelve:
#   $Phi_hat   : matriz [N x N*p] de coeficientes (sin intercepto)
#   $Sigma_hat : matriz de covarianza de residuos [N x N]
#   $lambdas   : lambda optimo por ecuacion (via cv.glmnet)
#   $sparsity  : fraccion de coeficientes llevados a cero
# -----------------------------------------------------------------------------
estimar_alasso_var <- function(Y, p = 1) {
  
  N <- ncol(Y)
  T <- nrow(Y)
  
  # Construir matrices de regresor y dependiente
  # Y_dep: T-p x N (variables dependientes)
  # X_reg: T-p x N*p (rezagos apilados) + columna de unos (intercepto)
  Y_dep <- Y[(p + 1):T, , drop = FALSE]
  
  X_lag_list <- lapply(1:p, function(lag) Y[(p + 1 - lag):(T - lag), , drop = FALSE])
  X_reg      <- do.call(cbind, X_lag_list)   # T-p x N*p (sin intercepto; glmnet lo agrega)
  
  T_eff <- nrow(Y_dep)
  
  # ------- Etapa 1: OLS para obtener pesos adaptativos -------
  # Agregar intercepto manualmente para OLS
  X_ols    <- cbind(1, X_reg)
  Phi_ols  <- solve(t(X_ols) %*% X_ols) %*% t(X_ols) %*% Y_dep  # (N*p+1) x N
  
  # Coeficientes OLS sin intercepto (fila 1 = intercepto, filas 2:(N*p+1) = Phi)
  coef_ols_sin_intercepto <- Phi_ols[2:(N * p + 1), , drop = FALSE]  # N*p x N
  
  # ------- Etapa 2: aLASSO ecuacion por ecuacion -------
  Phi_alasso <- matrix(0, nrow = N * p, ncol = N)
  residuos   <- matrix(0, nrow = T_eff, ncol = N)
  lambdas    <- numeric(N)
  
  for (j in 1:N) {
    y_j <- Y_dep[, j]
    
    # Pesos adaptativos: 1 / |phi_OLS_j| + epsilon para evitar division por cero
    # Un coeficiente OLS grande -> peso pequeno -> menos penalizacion -> se retiene
    # Un coeficiente OLS pequeno -> peso grande -> mas penalizacion -> tiende a cero
    eps         <- 1e-6
    w_adapt     <- 1 / (abs(coef_ols_sin_intercepto[, j]) + eps)
    
    # cv.glmnet con penalty.factor = pesos adaptativos
    # alpha = 1 -> LASSO puro (L1); con pesos adaptativos = adaptive LASSO
    set.seed(42)
    cv_fit <- cv.glmnet(
      x              = X_reg,
      y              = y_j,
      alpha          = 1,              # LASSO (L1)
      penalty.factor = w_adapt,        # pesos adaptativos por coeficiente
      intercept      = TRUE,
      standardize    = TRUE,           # estandarizar internamente (recomendado)
      nfolds         = 5               # 5-fold CV (adecuado para T pequeño)
    )
    
    lambda_opt   <- cv_fit$lambda.min  # lambda que minimiza MSE en CV
    lambdas[j]   <- lambda_opt
    
    coef_j       <- coef(cv_fit, s = lambda_opt)
    # coef_j[1] = intercepto, coef_j[2:(N*p+1)] = coeficientes de X_reg
    Phi_alasso[, j] <- as.numeric(coef_j)[2:(N * p + 1)]
    
    # Residuos de la ecuacion j
    residuos[, j] <- y_j - (X_reg %*% Phi_alasso[, j] + as.numeric(coef_j)[1])
  }
  
  # Matriz de covarianza de residuos (sin correccion de grados de libertad
  # porque con LASSO el numero efectivo de parametros es incierto)
  Sigma_hat <- (t(residuos) %*% residuos) / T_eff
  
  # Sparsity: fraccion de coeficientes en cero
  n_cero   <- sum(abs(Phi_alasso) < 1e-10)
  sparsity <- n_cero / (N * p * N)
  
  list(
    Phi_hat   = t(Phi_alasso),   # N x N*p (convencion: ecuacion en filas)
    Sigma_hat = Sigma_hat,
    lambdas   = lambdas,
    sparsity  = sparsity,
    residuos  = residuos
  )
}


# -----------------------------------------------------------------------------
# fevd_generalizada():
# Calcula la FEVD generalizada de Koop-Pesaran (DY 2012) a partir de
# Phi_hat (N x N*p) y Sigma (N x N).
# Devuelve theta_tilde: matriz N x N normalizada (filas suman 1).
# -----------------------------------------------------------------------------
fevd_generalizada <- function(Phi_hat, Sigma, H = 8, p = 1) {
  
  N <- nrow(Sigma)
  
  # Convertir VAR(p) companion form para obtener MA coefficients
  # Para VAR(1) simplificado: Psi_h = Phi^h
  # Para VAR(p) general se construye la companion matrix
  
  if (p == 1) {
    # Phi_hat es N x N directamente para p=1
    A <- Phi_hat   # N x N
  } else {
    # Companion matrix (N*p x N*p)
    A_top <- Phi_hat                                      # N x N*p
    A_bot <- cbind(diag(N * (p - 1)), matrix(0, N * (p - 1), N))
    A     <- rbind(A_top, A_bot)
  }
  
  # Calcular coeficientes MA: Psi_h = A^h (primeras N filas/columnas)
  Psi <- vector("list", H)
  A_h <- diag(nrow(A))   # A^0 = I
  for (h in 1:H) {
    A_h    <- A_h %*% A
    Psi[[h]] <- A_h[1:N, 1:N, drop = FALSE]
  }
  
  # FEVD generalizada (ecuacion 4 de DY 2012)
  sigma_diag <- diag(Sigma)   # varianzas propias (diagonal de Sigma)
  
  theta <- matrix(0, N, N)   # theta[i,j] = contribucion de j a varianza de i
  
  for (i in 1:N) {
    e_i <- matrix(0, N, 1); e_i[i] <- 1
    denom_i <- 0
    for (h in 0:(H - 1)) {
      Psi_h <- if (h == 0) diag(N) else Psi[[h]]
      denom_i <- denom_i + as.numeric(t(e_i) %*% Psi_h %*% Sigma %*% t(Psi_h) %*% e_i)
    }
    for (j in 1:N) {
      e_j    <- matrix(0, N, 1); e_j[j] <- 1
      numer  <- 0
      for (h in 0:(H - 1)) {
        Psi_h  <- if (h == 0) diag(N) else Psi[[h]]
        numer  <- numer + (as.numeric(t(e_i) %*% Psi_h %*% Sigma %*% e_j))^2
      }
      theta[i, j] <- numer / (sigma_diag[j] * denom_i)
    }
  }
  
  # Normalizar filas para que sumen 1
  theta_tilde <- theta / rowSums(theta)
  rownames(theta_tilde) <- colnames(theta_tilde) <- colnames(Phi_hat)
  theta_tilde
}


# -----------------------------------------------------------------------------
# indices_dy():
# Calcula TCI, TO, FROM, NET y la matriz pairwise a partir de theta_tilde.
# -----------------------------------------------------------------------------
indices_dy <- function(theta_tilde) {
  N   <- nrow(theta_tilde)
  
  # TO[i] = suma de contribuciones de i a otros / N  (excluye diagonal)
  #TO  <- colSums(theta_tilde - diag(diag(theta_tilde))) / N * 100
  
  # FROM[i] = suma de contribuciones que i recibe de otros / N
  #FROM <- rowSums(theta_tilde - diag(diag(theta_tilde))) / N * 100
  
  # NET[i] = TO[i] - FROM[i]
  #NET <- TO - FROM
  
  # TCI = media del TO (o equivalentemente, media del FROM)
  #TCI <- mean(TO)
  
  # version alineada con ConnectednessApproach
  TO   <- colSums(theta_tilde - diag(diag(theta_tilde))) * 100
  FROM <- rowSums(theta_tilde - diag(diag(theta_tilde))) * 100
  TCI  <- sum(theta_tilde - diag(diag(theta_tilde))) / N * 100
  
  # NET[i] = TO[i] - FROM[i]
  NET <- TO - FROM
  
  # Pairwise: theta_tilde tal cual (ya normalizada)
  list(TCI = TCI, TO = TO, FROM = FROM, NET = NET, pairwise = theta_tilde)
}


# =============================================================================
# PASO 4 - ADAPTIVE LASSO-VAR ESTATICO (MUESTRA COMPLETA)
# =============================================================================
cat("\n--- PASO 4: ADAPTIVE LASSO-VAR ESTATICO ---\n")

# Rezago seleccionado por BIC (mismo que el VAR OLS)
Y_data   <- coredata(Y_zoo)
sel_lags <- VARselect(Y_data, lag.max = 4, type = "const")
P_sel    <- sel_lags$selection["SC(n)"]
cat(sprintf("  Rezago seleccionado por BIC: p = %d\n", P_sel))

# Estimar aLASSO sobre muestra completa
alasso_full <- estimar_alasso_var(Y_data, p = P_sel)

cat(sprintf("  Sparsity del modelo: %.1f%% de coeficientes en cero\n",
            alasso_full$sparsity * 100))
cat("  Lambdas optimos por ecuacion:\n")
cat(paste0("    ", colnames(Y_data), ": ", round(alasso_full$lambdas, 5),
           collapse = "\n"), "\n")

# FEVD y indices estaticos
theta_full  <- fevd_generalizada(alasso_full$Phi_hat, alasso_full$Sigma_hat,
                                 H = 8, p = P_sel)
idx_full    <- indices_dy(theta_full)

# Tabla de connectedness estatica (formato similar a ConnectednessApproach$TABLE)
tabla_conn <- round(rbind(
  theta_full * 100,
  TO   = idx_full$TO,
  FROM = idx_full$FROM,
  NET  = idx_full$NET
), 2)

write.csv(tabla_conn,
          file.path(ruta_output, "lasso_bcch_04_connectedness_completa.csv"),
          row.names = TRUE)

cat("\n  Tabla de connectedness estatica (aLASSO):\n")
print(tabla_conn)


# =============================================================================
# PASO 5 - ROLLING-WINDOW ADAPTIVE LASSO-VAR
# =============================================================================
cat("\n--- PASO 5: ROLLING-WINDOW ADAPTIVE LASSO-VAR ---\n")

ventanas <- c(34)

resultados_tci_rw      <- list()
resultados_net_rw      <- list()
resultados_to_rw       <- list()
resultados_from_rw     <- list()
resultados_pairwise_rw <- list()
resultados_pw_rm_bancos <- list()
resultados_pw_rm_fp     <- list()

nombres_vec <- colnames(Y_data)   # vector de nombres de sectores

for (w in ventanas) {
  cat(sprintf("\n  Estimando rolling-window aLASSO (W = %d)...\n", w))
  
  T_total <- nrow(Y_data)
  n_ventanas <- T_total - w + 1
  
  # Preallocar listas de resultados
  tci_vec  <- numeric(n_ventanas)
  to_mat   <- matrix(NA, n_ventanas, length(nombres_vec))
  from_mat <- matrix(NA, n_ventanas, length(nombres_vec))
  net_mat  <- matrix(NA, n_ventanas, length(nombres_vec))
  
  # Indices para pairwise
  idx_bancos     <- which(nombres_vec == "Bancos")
  idx_fp         <- which(nombres_vec == "FP")
  idx_restomundo <- which(nombres_vec == "RestoMundo")
  
  pw_fp_bancos    <- numeric(n_ventanas)
  pw_bancos_fp    <- numeric(n_ventanas)
  pw_rm_bancos    <- numeric(n_ventanas)
  pw_bancos_rm    <- numeric(n_ventanas)
  pw_rm_fp        <- numeric(n_ventanas)
  pw_fp_rm        <- numeric(n_ventanas)
  
  for (t in 1:n_ventanas) {
    
    # Mostrar progreso cada 10 ventanas
    if (t %% 10 == 0) cat(sprintf("    Ventana %d / %d\r", t, n_ventanas))
    
    # Submuestra de la ventana
    Y_win <- Y_data[t:(t + w - 1), , drop = FALSE]
    
    # Estimar aLASSO en esta ventana
    # Capturar errores para ventanas problematicas sin detener el loop
    alasso_win <- tryCatch(
      estimar_alasso_var(Y_win, p = P_sel),
      error = function(e) NULL
    )
    
    if (is.null(alasso_win)) {
      # Si falla, propagar NA para esta ventana
      tci_vec[t]    <- NA
      to_mat[t, ]   <- NA
      from_mat[t, ] <- NA
      net_mat[t, ]  <- NA
      next
    }
    
    # FEVD y indices
    theta_win <- tryCatch(
      fevd_generalizada(alasso_win$Phi_hat, alasso_win$Sigma_hat,
                        H = 8, p = P_sel),
      error = function(e) NULL
    )
    
    if (is.null(theta_win)) {
      tci_vec[t]    <- NA
      to_mat[t, ]   <- NA
      from_mat[t, ] <- NA
      net_mat[t, ]  <- NA
      next
    }
    
    idx_win <- indices_dy(theta_win)
    
    tci_vec[t]    <- idx_win$TCI
    to_mat[t, ]   <- idx_win$TO
    from_mat[t, ] <- idx_win$FROM
    net_mat[t, ]  <- idx_win$NET
    
    # Pairwise: theta_win[i, j] = fraccion de varianza de i explicada por j
    # FP contagia a Bancos: cuanto de Bancos se explica por FP
    pw_fp_bancos[t] <- theta_win[idx_bancos, idx_fp]    * 100
    pw_bancos_fp[t] <- theta_win[idx_fp,     idx_bancos] * 100
    pw_rm_bancos[t] <- theta_win[idx_bancos, idx_restomundo] * 100
    pw_bancos_rm[t] <- theta_win[idx_restomundo, idx_bancos] * 100
    pw_rm_fp[t]     <- theta_win[idx_fp,     idx_restomundo] * 100
    pw_fp_rm[t]     <- theta_win[idx_restomundo, idx_fp]     * 100
  }
  
  # Fechas correspondientes a cada ventana (ultima fecha de cada ventana)
  fechas_rw <- index(Y_zoo)[w:T_total]
  
  cat(sprintf("\n    W=%d completado. TCI medio: %.2f%%\n", w, mean(tci_vec, na.rm = TRUE)))
  
  # --- Armar dataframes ---
  w_label <- paste0("W=", w)
  
  # TCI
  resultados_tci_rw[[w_label]] <- data.frame(
    fecha   = as.Date(fechas_rw),
    tci     = tci_vec,
    ventana = w_label
  )
  
  # NET
  colnames(net_mat) <- nombres_vec
  resultados_net_rw[[w_label]] <- as.data.frame(net_mat) %>%
    mutate(fecha = as.Date(fechas_rw), ventana = w_label) %>%
    pivot_longer(cols = -c(fecha, ventana), names_to = "sector", values_to = "value")
  
  # TO
  colnames(to_mat) <- nombres_vec
  resultados_to_rw[[w_label]] <- as.data.frame(to_mat) %>%
    mutate(fecha = as.Date(fechas_rw), ventana = w_label) %>%
    pivot_longer(cols = -c(fecha, ventana), names_to = "sector", values_to = "value")
  
  # FROM
  colnames(from_mat) <- nombres_vec
  resultados_from_rw[[w_label]] <- as.data.frame(from_mat) %>%
    mutate(fecha = as.Date(fechas_rw), ventana = w_label) %>%
    pivot_longer(cols = -c(fecha, ventana), names_to = "sector", values_to = "value")
  
  # Pairwise FP vs Bancos
  resultados_pairwise_rw[[w_label]] <- data.frame(
    fecha           = as.Date(fechas_rw),
    ventana         = w_label,
    `FP contagia a Bancos`    = pw_fp_bancos,
    `Bancos contagia a FP`    = pw_bancos_fp,
    check.names = FALSE
  ) %>%
    pivot_longer(
      cols      = c(`FP contagia a Bancos`, `Bancos contagia a FP`),
      names_to  = "direccion",
      values_to = "pairwise_value"
    )
  
  # Pairwise RestoMundo vs Bancos
  resultados_pw_rm_bancos[[w_label]] <- data.frame(
    fecha           = as.Date(fechas_rw),
    ventana         = w_label,
    `Resto Mundo contagia a Bancos` = pw_rm_bancos,
    `Bancos contagia a Resto Mundo` = pw_bancos_rm,
    check.names = FALSE
  ) %>%
    pivot_longer(
      cols      = c(`Resto Mundo contagia a Bancos`, `Bancos contagia a Resto Mundo`),
      names_to  = "direccion",
      values_to = "pairwise_value"
    )
  
  # Pairwise RestoMundo vs FP
  resultados_pw_rm_fp[[w_label]] <- data.frame(
    fecha           = as.Date(fechas_rw),
    ventana         = w_label,
    `Resto Mundo contagia a FP` = pw_rm_fp,
    `FP contagia a Resto Mundo` = pw_fp_rm,
    check.names = FALSE
  ) %>%
    pivot_longer(
      cols      = c(`Resto Mundo contagia a FP`, `FP contagia a Resto Mundo`),
      names_to  = "direccion",
      values_to = "pairwise_value"
    )
}

# --- Combinar y guardar CSVs ---
df_tci_rw       <- bind_rows(resultados_tci_rw)
df_net_rw       <- bind_rows(resultados_net_rw)
df_to_rw        <- bind_rows(resultados_to_rw)
df_from_rw      <- bind_rows(resultados_from_rw)
df_pairwise_rw  <- bind_rows(resultados_pairwise_rw)
df_pw_rm_bancos <- bind_rows(resultados_pw_rm_bancos)
df_pw_rm_fp     <- bind_rows(resultados_pw_rm_fp)

write_csv(df_tci_rw,       file.path(ruta_output, "lasso_bcch_05_tci_rolling_window.csv"))
write_csv(df_net_rw,       file.path(ruta_output, "lasso_bcch_06_net_rolling_window.csv"))
write_csv(df_to_rw,        file.path(ruta_output, "lasso_bcch_07_to_rolling_window.csv"))
write_csv(df_from_rw,      file.path(ruta_output, "lasso_bcch_08_from_rolling_window.csv"))
write_csv(df_pairwise_rw,  file.path(ruta_output, "lasso_bcch_09_pairwise_fp_bancos.csv"))
write_csv(df_pw_rm_bancos, file.path(ruta_output, "lasso_bcch_10_pairwise_RM_bancos.csv"))
write_csv(df_pw_rm_fp,     file.path(ruta_output, "lasso_bcch_11_pairwise_RM_FP.csv"))

cat("\n  CSVs guardados en salidas_lasso/\n")


# =============================================================================
# GRAFICOS (identicos en estructura al script VAR original)
# =============================================================================
cat("\n--- GRAFICOS ---\n")

# Paleta y tema comun
colores_ventanas <- c("W=32" = "gray50", "W=34" = "steelblue", "W=36" = "forestgreen")

tema_graficos <- theme_minimal(base_size = 10) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "gray95", color = NA),
    axis.text.x      = element_text(angle = 40, hjust = 1)
  )

# Orden de narración
orden_narrativo <- c("FP", "Bancos", "BancoCentral", "Seguros", 
                     "FMM", "FNM", "OFIs", "RestoMundo")
#orden_original <- unname(nombres_sectores[orden_sectores])

df_net_rw <- df_net_rw %>% mutate(sector = factor(sector, levels = orden_narrativo))
df_to_rw  <- df_to_rw %>% mutate(sector = factor(sector, levels = orden_narrativo))
df_from_rw  <- df_from_rw %>% mutate(sector = factor(sector, levels = orden_narrativo))



# Funcion auxiliar para agregar eventos comunes a todos los graficos
agregar_eventos <- function(p) {
  p +
    geom_vline(xintercept = as.Date("2008-09-01"), color = "gray40", linetype = "dotted") +
    geom_vline(xintercept = as.Date("2019-10-01"), color = "gray40", linetype = "dotted") +
    #retiros
    geom_vline(xintercept = as.Date("2020-07-01"), color = "red",
               linetype = "dotted", linewidth = 0.4) +
    geom_vline(xintercept = as.Date("2020-12-01"), color = "red",
               linetype = "dotted", linewidth = 0.4) +
    geom_vline(xintercept = as.Date("2021-04-01"), color = "red",
               linetype = "dotted", linewidth = 0.4) +
    annotate("rect",
             xmin = as.Date("2020-07-01"), xmax = as.Date("2021-12-31"),
             ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red") +
    annotate("text", x = as.Date("2021-04-01"), y = Inf, label = "Periodo\nRetiros", 
             vjust = 1.5, hjust = 0.5, size = 3, color = "darkred", fontface = "bold") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y-T%q") +
    scale_color_manual(values = colores_ventanas) +
    scale_linetype_manual(values = c("W=32" = "dashed", "W=34" = "solid", "W=36" = "dotdash"))
}

# --- Grafico 1: TCI ---
p_tci <- ggplot(df_tci_rw,
                aes(x = fecha, y = tci, color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.8) +
  labs(
    title    = "Total Connectivity Index (TCI) — Adaptive LASSO-VAR",
    subtitle = "Fraccion de la varianza de error de prediccion explicada por shocks cruzados",
    x        = "Periodo", y = "TCI (%)",
    color    = "Ventana", linetype = "Ventana"
  ) +
  tema_graficos
p_tci <- agregar_eventos(p_tci)

# --- Grafico 2: NET ---
p_net <- ggplot(df_net_rw,
                aes(x = fecha, y = value, color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ sector, scales = "free_y", ncol = 2) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  labs(
    title    = "Net Directional Connectedness por Sector — Adaptive LASSO-VAR",
    subtitle = "Valores > 0: Transmisor Neto | Valores < 0: Receptor Neto",
    x        = "Periodo", y = "Net Connectedness",
    color    = "Ventana", linetype = "Ventana"
  ) +
  tema_graficos
p_net <- agregar_eventos(p_net)

# --- Grafico 3: TO ---
p_to <- ggplot(df_to_rw,
               aes(x = fecha, y = value, color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ sector, scales = "free_y", ncol = 2) +
  labs(
    title    = "To-Connectedness por Sector — Adaptive LASSO-VAR",
    subtitle = "Proporcion de la varianza de otros sectores explicada por shocks de este sector",
    x        = "Periodo", y = "To-Connectedness (%)",
    color    = "Ventana", linetype = "Ventana"
  ) +
  tema_graficos
p_to <- agregar_eventos(p_to)

# --- Grafico 4: FROM ---
p_from <- ggplot(df_from_rw,
                 aes(x = fecha, y = value, color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ sector, scales = "free_y", ncol = 2) +
  labs(
    title    = "From-Connectedness por Sector — Adaptive LASSO-VAR",
    subtitle = "Proporcion de la varianza del sector explicada por shocks de otros sectores",
    x        = "Periodo", y = "From-Connectedness (%)",
    color    = "Ventana", linetype = "Ventana"
  ) +
  tema_graficos
p_from <- agregar_eventos(p_from)

# --- Grafico 5: Pairwise FP vs Bancos ---
p_pairwise <- ggplot(df_pairwise_rw,
                     aes(x = fecha, y = pairwise_value,
                         color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ direccion, scales = "free_y", ncol = 1) +
  labs(
    title    = "Conectividad Bilateral: Bancos vs Fondos de Pensiones — Adaptive LASSO-VAR",
    subtitle = "Varianza explicada (%) en cada direccion",
    x        = "Periodo", y = "Varianza Explicada (%)",
    color    = "Ventana", linetype = "Ventana"
  ) +
  tema_graficos
p_pairwise <- agregar_eventos(p_pairwise)

# --- Grafico 6: Pairwise RestoMundo vs Bancos ---
p_pw_rm_bancos <- ggplot(df_pw_rm_bancos,
                         aes(x = fecha, y = pairwise_value,
                             color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ direccion, scales = "free_y", ncol = 1) +
  labs(
    title    = "Conectividad Bilateral: Resto Mundo vs Bancos — Adaptive LASSO-VAR",
    subtitle = "Varianza explicada (%) en cada direccion",
    x        = "Periodo", y = "Varianza Explicada (%)",
    color    = "Ventana", linetype = "Ventana"
  ) +
  tema_graficos
p_pw_rm_bancos <- agregar_eventos(p_pw_rm_bancos)

# --- Grafico 7: Pairwise RestoMundo vs FP ---
p_pw_rm_fp <- ggplot(df_pw_rm_fp,
                     aes(x = fecha, y = pairwise_value,
                         color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ direccion, scales = "free_y", ncol = 1) +
  labs(
    title    = "Conectividad Bilateral: Resto Mundo vs FP — Adaptive LASSO-VAR",
    subtitle = "Varianza explicada (%) en cada direccion",
    x        = "Periodo", y = "Varianza Explicada (%)",
    color    = "Ventana", linetype = "Ventana"
  ) +
  tema_graficos
p_pw_rm_fp <- agregar_eventos(p_pw_rm_fp)

# --- Guardar graficos ---
ggsave(file.path(ruta_output, "figs", "lasso_bcch_01_tci_rolling_window.png"),
       p_tci,         width = 10, height = 6,  dpi = 150, bg = "white")
ggsave(file.path(ruta_output, "figs", "lasso_bcch_02_net_rolling_window.png"),
       p_net,         width = 10, height = 8,  dpi = 150, bg = "white")
ggsave(file.path(ruta_output, "figs", "lasso_bcch_03_to_rolling_window.png"),
       p_to,          width = 10, height = 8,  dpi = 150, bg = "white")
ggsave(file.path(ruta_output, "figs", "lasso_bcch_04_from_rolling_window.png"),
       p_from,        width = 10, height = 8,  dpi = 150, bg = "white")
ggsave(file.path(ruta_output, "figs", "lasso_bcch_05_pairwise_fp_bancos.png"),
       p_pairwise,    width = 10, height = 8,  dpi = 150, bg = "white")
ggsave(file.path(ruta_output, "figs", "lasso_bcch_06_pairwise_rm_bancos.png"),
       p_pw_rm_bancos, width = 10, height = 8, dpi = 150, bg = "white")
ggsave(file.path(ruta_output, "figs", "lasso_bcch_07_pairwise_RM_FP.png"),
       p_pw_rm_fp,    width = 10, height = 8,  dpi = 150, bg = "white")

cat("\n  Graficos 1 a 7 guardados en salidas_lasso_bcch/figs/\n")


# =============================================================================
# TABLA RESUMEN 4 FECHAS (misma logica que el script VAR original)
# =============================================================================
cat("\n--- TABLA RESUMEN 4 FECHAS ---\n")

fechas_target <- as.Date(c("2019-10-01", "2020-04-01", "2020-10-01", "2022-01-01"))
periodos_lbl  <- c("2019T4", "2020T2", "2020T4", "2022T1")

# Usar W=34 como ventana principal (misma convencion que el script VAR)
tci_tab <- df_tci_rw %>%
  filter(ventana == "W=34", fecha %in% fechas_target) %>%
  select(fecha, TCI_W34 = tci)

from_fp_tab <- df_from_rw %>%
  filter(ventana == "W=34", sector == "FP", fecha %in% fechas_target) %>%
  select(fecha, From_FP_W34 = value)

to_fp_tab <- df_to_rw %>%
  filter(ventana == "W=34", sector == "FP", fecha %in% fechas_target) %>%
  select(fecha, To_FP_W34 = value)

pw_fp_bancos_tab <- df_pairwise_rw %>%
  filter(ventana == "W=34", direccion == "FP contagia a Bancos",
         fecha %in% fechas_target) %>%
  select(fecha, Pairwise_FP_Bancos = pairwise_value)

pw_fp_row_tab <- df_pw_rm_fp %>%
  filter(ventana == "W=34", direccion == "FP contagia a Resto Mundo",
         fecha %in% fechas_target) %>%
  select(fecha, Pairwise_FP_RoW = pairwise_value)

pw_bancos_fp_tab <- df_pairwise_rw %>%
  filter(ventana == "W=34", direccion == "Bancos contagia a FP",
         fecha %in% fechas_target) %>%
  select(fecha, Pairwise_Bancos_FP = pairwise_value)

tabla_resumen <- list(
  tci_tab, from_fp_tab, to_fp_tab,
  pw_fp_bancos_tab, pw_fp_row_tab, pw_bancos_fp_tab
) %>%
  reduce(left_join, by = "fecha") %>%
  mutate(
    periodo    = periodos_lbl[match(fecha, fechas_target)],
    Net_FP_W34 = To_FP_W34 - From_FP_W34
  ) %>%
  select(periodo, fecha,
         TCI_W34, From_FP_W34, To_FP_W34, Net_FP_W34,
         Pairwise_FP_Bancos, Pairwise_FP_RoW, Pairwise_Bancos_FP)

print(tabla_resumen)

write_csv(tabla_resumen,  file.path(ruta_output, "lasso_bcch_TABLA_RESUMEN_4FECHAS.csv"))
write_xlsx(tabla_resumen, file.path(ruta_output, "lasso_bcch_TABLA_RESUMEN_4FECHAS.xlsx"))

cat("\n===============================================================\n")
cat("PROCESO COMPLETADO. OUTPUTS EN: salidas_lasso/\n")
cat("===============================================================\n")


# =============================================================================
# NOTA: ALTERNATIVA SIMPLE CON ConnectednessApproach model = "LASSO"
# =============================================================================
# Si prefieres usar directamente el LASSO estandar del paquete (sin adaptive
# weights), puedes reemplazar el loop de ventanas por esto:
#
#   dy_rw_lasso <- ConnectednessApproach(
#     x             = Y_zoo,
#     nlag          = P_sel,
#     nfore         = 8,
#     model         = "LASSO",      # LASSO estandar via BigVAR
#     connectedness = "Time",
#     window.size   = 32
#   )
#
# Ventaja: mas simple, mismo output que VAR irrestricto, sin codigo adicional.
# Desventaja: LASSO estandar introduce sesgo en coeficientes grandes.
#             El adaptive LASSO implementado arriba corrige ese sesgo.
# =============================================================================