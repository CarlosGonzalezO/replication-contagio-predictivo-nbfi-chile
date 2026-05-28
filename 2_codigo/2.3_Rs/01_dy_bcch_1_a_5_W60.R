# =============================================================================
# TESIS MAGISTER ECONOMIA UC
# Intermediacion financiera no bancaria e interconexiones sectoriales en Chile
# =============================================================================
# Script: 01_dy_bcch_1_a_5_W60.R
# Autor:  Carlos Gonzalez
# Fecha:  Mayo 2026
#
# OBJETIVO:
# Implementar los Pasos 1 a 5 del analisis Diebold-Yilmaz:
#   Paso 1: Setup del entorno (instalacion de paquetes, carga de librerias)
#   Paso 2: Carga y validacion de la data (panel_var_dy.csv)
#   Paso 3: Tests de estacionariedad (ADF y KPSS por sector)
#   Paso 4: VAR estatico sobre muestra completa + matriz de connectedness
#   Paso 5: Rolling-window complementario (TCI y Sectorial, ventanas 20, 24, 28)
#
# INPUTS:
#   panel_var_dy.csv (output del pipeline Python consolidado_cnsi.py)
#
# OUTPUTS GENERADOS:
#   /salidas/01_serie_var_wide.csv         data en formato wide para VAR
#   /salidas/02_tests_estacionariedad.csv  resultados ADF y KPSS
#   /salidas/03_var_completo_resumen.txt   resumen del VAR estimado
#   /salidas/04_connectedness_completa.csv tabla de connectedness DY
#   /salidas/05_tci_rolling_window.csv     datos del TCI dinamico (20, 24, 28)
#   /salidas/06_net_rolling_window.csv     datos Net Connectedness por sector
#   /salidas/figs/                         graficos de validacion y rolling
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
ruta_output <- ruta / "3.2_resultado_dy/salidas_W60"

# =============================================================================
# PASO 1 - SETUP DEL ENTORNO
# =============================================================================
paquetes <- c(
  "ConnectednessApproach", "vars", "urca", "zoo", 
  "tidyverse", "lubridate", "knitr", "kableExtra"
)
paquetes_faltantes <- paquetes[!paquetes %in% installed.packages()[, "Package"]]
if (length(paquetes_faltantes) > 0) install.packages(paquetes_faltantes)

library(ConnectednessApproach)
library(vars)
library(urca)
library(zoo)
library(tidyverse)
library(lubridate)
library(tseries)

dir.create(ruta_output,           showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ruta_output, "figs"), showWarnings = FALSE)

options(scipen = 999)
set.seed(42)

cat("\n===============================================================\n")
cat("DIEBOLD-YILMAZ: PASOS 1 A 5\n")
cat("===============================================================\n")


# =============================================================================
# PASO 2 - CARGA Y VALIDACION DE DATA
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

write_csv(dy_wide, file.path(ruta_output, "01_serie_var_wide.csv"))

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
  year <- as.numeric(substr(periodo, 1, 4))
  trim <- as.numeric(substr(periodo, 6, 6))
  mes_inicio <- (trim - 1) * 3 + 1
  as.Date(paste0(year, "-", sprintf("%02d", mes_inicio), "-01"))
}

fechas <- periodo_a_fecha(dy_wide$periodo)
Y_matrix <- as.matrix(dy_wide[, orden_sectores])
colnames(Y_matrix) <- nombres_sectores[orden_sectores]
Y_zoo <- zoo(Y_matrix, order.by = fechas)

# -------------------- Validacion visual exploratoria --------------------
# Plot de las series para detectar outliers y patrones antes del VAR.

dy_largo_para_plot <- dy_wide %>%
  select(periodo, all_of(orden_sectores)) %>%   # <-- FIX: solo sectores del VAR
  pivot_longer(-periodo, names_to = "sector_codigo", values_to = "d_log_A") %>%
  mutate(
    sector_nombre = nombres_sectores[sector_codigo],
    sector_nombre = factor(sector_nombre,        # <-- ordena los paneles
                           levels = unname(nombres_sectores[orden_sectores])),
    fecha = as.Date(paste0(substr(periodo, 1, 4), "-",
                           (as.numeric(substr(periodo, 6, 6)) - 1) * 3 + 1, "-01"))
  )

p_series <- ggplot(dy_largo_para_plot,
                   aes(x = fecha, y = d_log_A, color = sector_nombre)) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~ sector_nombre, scales = "free_y", ncol = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  # Bandas verticales para los tres episodios.
  geom_vline(xintercept = as.Date("2008-09-01"), color = "gray40",
             linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = as.Date("2019-10-01"), color = "gray40",
             linetype = "dotted", linewidth = 0.4) +
  #retiros
  geom_vline(xintercept = as.Date("2020-07-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = as.Date("2020-12-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = as.Date("2021-04-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  annotate("rect", xmin = as.Date("2020-07-01"), xmax = as.Date("2021-12-31"),
           ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red3") +
  labs(
    title = expression(paste(Delta, " log ", "Real Notional Stock",
                             " - serie del VAR Diebold-Yilmaz")),
    subtitle = "Lineas: GFC 2008T3, Estallido 2019T4 | Banda roja: Retiros 2020T3-2021T4",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "gray95", color = NA))

ggsave(file.path(ruta_output, "figs", "bcch_00_series_input_var.png"),
       p_series, width = 10, height = 6, dpi = 150)
cat("Figura guardada: bcch_00_series_input_var.png\n")


# =============================================================================
# PASO 3 - TESTS DE ESTACIONARIEDAD
# =============================================================================
cat("\n--- PASO 3: TESTS DE ESTACIONARIEDAD (ADF, PP, KPSS y Zivot-Andrews) ---\n")

# Helper: convierte la posicion del break (entero) a etiqueta trimestral YYYYTq.
# La serie del VAR comienza en 2003T2 (filter periodo != "2003T1" en Paso 2),
# por tanto la posicion 1 corresponde a 2003T2.
posicion_a_trimestre <- function(pos, inicio_anio = 2003, inicio_trim = 2) {
  if (is.na(pos) || pos < 1) return(NA_character_)
  trim_total <- (inicio_anio * 4 + (inicio_trim - 1)) + (pos - 1)
  anio_out   <- trim_total %/% 4
  trim_out   <- (trim_total %% 4) + 1
  sprintf("%dT%d", anio_out, trim_out)
}

testear_estacionariedad <- function(serie, nombre_serie) {
  serie_limpia <- na.omit(serie)
  
  # --- ADF con drift (constante); H0: raiz unitaria ---
  adf <- ur.df(serie_limpia, type = "drift", selectlags = "AIC")
  adf_stat   <- as.numeric(adf@teststat[1])
  adf_cv5    <- as.numeric(adf@cval[1, "5pct"])
  adf_rechaza_H0 <- adf_stat < adf_cv5   # rechazar = ESTACIONARIA
  
  # --- Phillips-Perron con drift (constante); H0: raiz unitaria ---
  # Test no parametrico, robusto a autocorrelacion y heterocedasticidad de los residuos.
  pp <- ur.pp(serie_limpia, type = "Z-tau", model = "constant", lags = "short")
  pp_stat <- as.numeric(pp@teststat)
  pp_cv5  <- as.numeric(pp@cval[1, "5pct"])
  pp_rechaza_H0 <- pp_stat < pp_cv5      # rechazar = ESTACIONARIA
  
  # --- KPSS con regresion sobre constante (level); H0: estacionaria ---
  kpss_l <- ur.kpss(serie_limpia, type = "mu", lags = "short")
  kpss_l_stat <- as.numeric(kpss_l@teststat)
  kpss_l_cv5  <- as.numeric(kpss_l@cval[1, "5pct"])
  kpss_l_rechaza_H0 <- kpss_l_stat > kpss_l_cv5  # rechazar = NO ESTACIONARIA
  
  # --- Zivot-Andrews con break endogeno en intercepto y tendencia ("both") ---
  # H0: raiz unitaria con drift; H1: estacionaria con un quiebre estructural.
  # Permite distinguir entre raiz unitaria genuina y rechazo espurio de KPSS por breaks.
  za <- tryCatch(
    ur.za(serie_limpia, model = "both", lag = 1),
    error = function(e) NULL
  )
  if (!is.null(za)) {
    za_stat       <- as.numeric(za@teststat)
    # ur.za@cval es un vector nombrado de longitud 3 (no matriz).
    # Los nombres pueden variar entre versiones de urca: probar varios alias.
    cv_vec <- za@cval
    nm <- names(cv_vec)
    idx5 <- which(nm %in% c("5pct", "5%", "5% critical value"))
    if (length(idx5) == 1) {
      za_cv5 <- as.numeric(cv_vec[idx5])
    } else {
      # Fallback: el critico al 5% siempre es el segundo elemento (1%, 5%, 10%).
      za_cv5 <- as.numeric(cv_vec[2])
    }
    za_rechaza_H0 <- za_stat < za_cv5      # rechazar = ESTACIONARIA c/ break
    za_break_pos  <- as.integer(za@bpoint)
    za_break_lbl  <- posicion_a_trimestre(za_break_pos)
  } else {
    za_stat <- NA_real_; za_cv5 <- NA_real_
    za_rechaza_H0 <- NA; za_break_pos <- NA_integer_
    za_break_lbl <- NA_character_
  }
  
  # --- Conclusion informativa basada en los 4 tests ---
  # Logica: ADF + PP convergentes pesan mas que KPSS (sensible a breaks).
  # ZA es el desempate cuando KPSS rechaza pero ADF/PP no.
  ru_unitaria_consistente <- adf_rechaza_H0 & pp_rechaza_H0
  
  conclusion <- if (ru_unitaria_consistente & !kpss_l_rechaza_H0) {
    "ESTACIONARIA"
  } else if (ru_unitaria_consistente & kpss_l_rechaza_H0) {
    if (isTRUE(za_rechaza_H0)) "ESTACIONARIA c/ BREAK ESTRUCTURAL"
    else "REVISAR (conflicto KPSS sin break detectado)"
  } else if (!ru_unitaria_consistente & kpss_l_rechaza_H0) {
    "NO ESTACIONARIA"
  } else {
    "REVISAR"
  }
  
  data.frame(
    sector_codigo         = nombre_serie,
    sector_nombre         = nombres_sectores[nombre_serie],
    ADF_stat              = round(adf_stat, 3),
    ADF_cv5               = round(adf_cv5, 3),
    ADF_rechaza_H0        = adf_rechaza_H0,
    PP_stat               = round(pp_stat, 3),
    PP_cv5                = round(pp_cv5, 3),
    PP_rechaza_H0         = pp_rechaza_H0,
    KPSS_lvl_stat         = round(kpss_l_stat, 3),
    KPSS_lvl_cv5          = round(kpss_l_cv5, 3),
    KPSS_lvl_rechaza_H0   = kpss_l_rechaza_H0,
    ZA_stat               = round(za_stat, 3),
    ZA_cv5                = round(za_cv5, 3),
    ZA_rechaza_H0         = za_rechaza_H0,
    ZA_break_trimestre    = za_break_lbl,
    conclusion            = conclusion,
    stringsAsFactors      = FALSE
  )
}

resultados_estacionariedad <- map_dfr(orden_sectores, function(s) {
  testear_estacionariedad(coredata(Y_zoo)[, nombres_sectores[s]], s)
})

# Imprimir en consola para revision inmediata
cat("\n")
print(resultados_estacionariedad %>% 
        select(sector_nombre, ADF_stat, PP_stat, KPSS_lvl_stat,
               ZA_stat, ZA_break_trimestre, conclusion))
cat(sprintf(
  "\nValores criticos al 5%%: ADF = %.3f | PP = %.3f | KPSS-level = %.3f | ZA(both) = %.3f\n",
  resultados_estacionariedad$ADF_cv5[1],
  resultados_estacionariedad$PP_cv5[1],
  resultados_estacionariedad$KPSS_lvl_cv5[1],
  resultados_estacionariedad$ZA_cv5[1]
))

write_csv(resultados_estacionariedad, file.path(ruta_output, "bcch_02_tests_estacionariedad.csv"))

# =============================================================================
# PASO 4 - VAR ESTATICO SOBRE MUESTRA COMPLETA + CONNECTEDNESS DY
# =============================================================================
cat("\n--- PASO 4: VAR + CONNECTEDNESS SOBRE MUESTRA COMPLETA ---\n")

Y_data_var <- coredata(Y_zoo)
seleccion_lags <- VARselect(Y_data_var, lag.max = 4, type = "const")
P_seleccionado <- seleccion_lags$selection["SC(n)"]

var_completo <- VAR(Y_data_var, p = P_seleccionado, type = "const")

sink(file.path(ruta_output, "03_var_completo_resumen.txt"))
print(summary(var_completo))
sink()

dy_resultado <- ConnectednessApproach(
  x             = Y_zoo,
  nlag          = P_seleccionado,
  nfore         = 8,
  model         = "VAR",
  connectedness = "Time",
  window.size   = NULL
)

write.csv(dy_resultado$TABLE, file.path(ruta_output, "04_connectedness_completa.csv"), row.names = TRUE)


# =============================================================================
# PASO 5 - ROLLING-WINDOW CONNECTEDNESS (TCI, NET, TO y FROM Sectorial)
# =============================================================================
cat("\n--- PASO 5: ROLLING-WINDOW (TCI, NET, TO y FROM) ---\n")

ventanas <-  c(60,62,64)
resultados_tci_rw  <- list()
resultados_net_rw  <- list()
resultados_to_rw   <- list()
resultados_from_rw <- list() # Nueva lista para FROM
resultados_pairwise_rw <- list() # Bancos FP
resultados_pairwise_rm_Bancos <- list() # RestoMundo Bancos 
resultados_pairwise_rm_FP <- list() # RestoMundo FP
resultados_pairwise_fp_bcch <- list() # FP BC

for (w in ventanas) {
  cat(sprintf("Estimando VAR rolling-window (W = %d)...\n", w))
  
  dy_rw <- ConnectednessApproach(
    x             = Y_zoo,
    nlag          = P_seleccionado,
    nfore         = 8,
    model         = "VAR",
    connectedness = "Time",
    window.size   = w
  )
  
  fechas_totales <- index(Y_zoo)
  fechas_rw <- tail(fechas_totales, length(as.numeric(dy_rw$TCI)))
  
  # 1. Extraer TCI
  resultados_tci_rw[[as.character(w)]] <- data.frame(
    fecha = as.Date(fechas_rw),
    tci = as.numeric(dy_rw$TCI),
    ventana = paste0("W=", w)
  )
  
  # 2. Extraer Net Connectedness
  net_matrix <- as.matrix(dy_rw$NET)
  colnames(net_matrix) <- unname(nombres_sectores[orden_sectores])
  resultados_net_rw[[as.character(w)]] <- as.data.frame(net_matrix) %>%
    mutate(fecha = as.Date(fechas_rw), ventana = paste0("W=", w)) %>%
    pivot_longer(cols = -c(fecha, ventana), names_to = "sector", values_to = "value")
  
  # 3. Extraer TO Connectedness (Hacia otros)
  to_matrix <- as.matrix(dy_rw$TO)
  colnames(to_matrix) <- unname(nombres_sectores[orden_sectores])
  resultados_to_rw[[as.character(w)]] <- as.data.frame(to_matrix) %>%
    mutate(fecha = as.Date(fechas_rw), ventana = paste0("W=", w)) %>%
    pivot_longer(cols = -c(fecha, ventana), names_to = "sector", values_to = "value")
  
  # 4. Extraer FROM Connectedness (Desde otros)
  from_matrix <- as.matrix(dy_rw$FROM)
  colnames(from_matrix) <- unname(nombres_sectores[orden_sectores])
  resultados_from_rw[[as.character(w)]] <- as.data.frame(from_matrix) %>%
    mutate(fecha = as.Date(fechas_rw), ventana = paste0("W=", w)) %>%
    pivot_longer(cols = -c(fecha, ventana), names_to = "sector", values_to = "value")
  
  # 5. Extraer Pairwise Connectedness (Bancos vs FP)
  # dy_rw$CT es un array 3D con dimensiones: [Receptor, Transmisor, Tiempo]
  idx_bancos <- which(unname(nombres_sectores[orden_sectores]) == "Bancos")
  idx_fp     <- which(unname(nombres_sectores[orden_sectores]) == "FP")
  
  # 6. Extraer Pairwise Connectedness (RestoMundo vs Banco)
  idx_bancos <- which(unname(nombres_sectores[orden_sectores]) == "Bancos")
  idx_RestoMundo <- which(unname(nombres_sectores[orden_sectores]) == "RestoMundo")
  
  # 7. Extraer Pairwise Connectedness (RestoMundo vs FP)
  idx_FP <- which(unname(nombres_sectores[orden_sectores]) == "FP")
  idx_RestoMundo <- which(unname(nombres_sectores[orden_sectores]) =="RestoMundo")
  
  # Cuanta varianza de Bancos es explicada por shocks de FP (FP contagia a Bancos)
  fp_to_bancos <- dy_rw$CT[idx_bancos, idx_fp, ]
  
  # Cuanta varianza de FP es explicada por shocks de Bancos (Bancos contagia a FP)
  bancos_to_fp <- dy_rw$CT[idx_fp, idx_bancos, ]
  
  # Cuanta varianza de Bancos es explicada por shocks de RestoMundo (RestoMundo contagia a Bancos)
  RestoMundo_to_Bancos <- dy_rw$CT[idx_bancos,idx_RestoMundo, ]
  
  # Cuanta varianza de RestoMundo es exlicada por shocks de Bancos (Bancos contagia a RestoMundo)
  Bancos_to_RestoMundo <- dy_rw$CT[idx_RestoMundo,idx_bancos, ]
  
  # Cuanta varianza de FP es explicada por shocks de RestoMundo (RestoMundo contagia a FP) 
  RestoMundo_to_FP <- dy_rw$CT[idx_FP,idx_RestoMundo , ]
  
  # Cuanta varianza de RestoMundo es explicada por shocks de FP (FP contagia a RestoMundo)
  FP_to_RestoMundo <- dy_rw$CT[idx_RestoMundo,idx_FP, ]
  
  resultados_pairwise_rw[[as.character(w)]] <- data.frame(
    fecha = as.Date(fechas_rw),
    ventana = paste0("W=", w),
    `FP contagia a Bancos` = as.numeric(fp_to_bancos) *100,
    `Bancos contagia a FP` = as.numeric(bancos_to_fp) *100,
    check.names = FALSE
  ) %>%
    pivot_longer(
      cols = c(`FP contagia a Bancos`, `Bancos contagia a FP`),
      names_to = "direccion",
      values_to = "pairwise_value"
    )
  
  resultados_pairwise_rm_Bancos[[as.character(w)]] <- data.frame(
    fecha = as.Date(fechas_rw),
    ventana = paste0("W=", w),
    `Resto Mundo contagia a Bancos` = as.numeric(RestoMundo_to_Bancos)*100,
    `Bancos contagia a Resto Mundo` = as.numeric(Bancos_to_RestoMundo)*100,
    check.names = FALSE
  ) %>%
    pivot_longer(
      cols = c(`Resto Mundo contagia a Bancos`, `Bancos contagia a Resto Mundo`),
      names_to = "direccion",
      values_to = "pairwise_value"
    )
  
  
  resultados_pairwise_rm_FP[[as.character(w)]] <- data.frame(
    fecha = as.Date(fechas_rw),
    ventana = paste0("W=", w),
    `Resto Mundo contagia a FP` = as.numeric(RestoMundo_to_FP)*100, # <-- Corregido
    `FP contagia a Resto Mundo` = as.numeric(FP_to_RestoMundo)*100,
    check.names = FALSE
  ) %>%
    pivot_longer(
      cols = c(`Resto Mundo contagia a FP`, `FP contagia a Resto Mundo`), # <-- Corregido
      names_to = "direccion",
      values_to = "pairwise_value"
    )
  # 8. Extraer Pairwise Connectedness (FP vs BancoCentral)
  idx_FP <- which(unname(nombres_sectores[orden_sectores]) == "FP")
  idx_BancoCentral <- which(unname(nombres_sectores[orden_sectores]) == "BancoCentral")
  
  # Cuanta varianza de BancoCentral es explicada por shocks de FP (FP contagia a BancoCentral)
  FP_to_BancoCentral <- dy_rw$CT[idx_BancoCentral, idx_FP, ]
  
  # Cuanta varianza de FP es explicada por shocks de BancoCentral (BancoCentral contagia a FP)
  BancoCentral_to_FP <- dy_rw$CT[idx_FP, idx_BancoCentral, ]
  
  resultados_pairwise_fp_bcch[[as.character(w)]] <- data.frame(
    fecha = as.Date(fechas_rw),
    ventana = paste0("W=", w),
    `FP contagia a BancoCentral` = as.numeric(FP_to_BancoCentral) * 100,
    `BancoCentral contagia a FP` = as.numeric(BancoCentral_to_FP) * 100,
    check.names = FALSE
  ) %>%
    pivot_longer(
      cols = c(`FP contagia a BancoCentral`, `BancoCentral contagia a FP`),
      names_to = "direccion",
      values_to = "pairwise_value"
    )
}

# Combinar y guardar CSVs
df_tci_rw  <- bind_rows(resultados_tci_rw)
df_net_rw  <- bind_rows(resultados_net_rw)
df_to_rw   <- bind_rows(resultados_to_rw)
df_from_rw <- bind_rows(resultados_from_rw)
df_pairwise_rw <- bind_rows(resultados_pairwise_rw)
df_pairwise_rm_Bancos <- bind_rows(resultados_pairwise_rm_Bancos)
df_pairwise_rm_FP <- bind_rows(resultados_pairwise_rm_FP)
df_pairwise_fp_bcch <- bind_rows(resultados_pairwise_fp_bcch)


write_csv(df_tci_rw,  file.path(ruta_output, "bcch_05_tci_rolling_window.csv"))
write_csv(df_net_rw,  file.path(ruta_output, "bcch_06_net_rolling_window.csv"))
write_csv(df_to_rw,   file.path(ruta_output, "bcch_07_to_rolling_window.csv"))
write_csv(df_from_rw, file.path(ruta_output, "bcch_08_from_rolling_window.csv"))
write_csv(df_pairwise_rw, file.path(ruta_output, "bcch_09_pairwise_fp_bancos.csv"))
write_csv(df_pairwise_rm_Bancos, file.path(ruta_output, "bcch_10_pairwise_RM_bancos.csv"))
write_csv(df_pairwise_rm_FP, file.path(ruta_output, "bcch_11_pairwise_RM_FP.csv"))
write_csv(df_pairwise_fp_bcch, file.path(ruta_output, "bcch_12_pairwise_fp_bcch.csv"))

# --- GRAFICOS ---
tema_graficos <- theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray95", color = NA),
    strip.text = element_text(face = "bold"), # Pone los nombres de los sectores en negrita
    axis.text.x = element_text(angle = 40, hjust = 1, color = "gray30"),
    axis.text.y = element_text(color = "gray30"),
    panel.grid.minor = element_blank(), # Elimina las líneas de cuadrícula menores para no saturar
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "gray40", size = 10)
  )

# Orden de narración
orden_narrativo <- c("FP", "Bancos", "BancoCentral", "Seguros", 
                     "FMM", "FNM", "OFIs", "RestoMundo")
#orden_original <- unname(nombres_sectores[orden_sectores])

df_net_rw <- df_net_rw %>% mutate(sector = factor(sector, levels = orden_narrativo))
df_to_rw  <- df_to_rw %>% mutate(sector = factor(sector, levels = orden_narrativo))
df_from_rw  <- df_from_rw %>% mutate(sector = factor(sector, levels = orden_narrativo))

# Grafico 1: TCI
p_tci <- ggplot(df_tci_rw, aes(x = fecha, y = tci, color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = as.Date("2008-09-01"), color = "gray40", linetype = "dotted") +
  geom_vline(xintercept = as.Date("2019-10-01"), color = "gray40", linetype = "dotted") +
  #retiros
  geom_vline(xintercept = as.Date("2020-07-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = as.Date("2020-12-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = as.Date("2021-04-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  annotate("rect", xmin = as.Date("2020-07-01"), xmax = as.Date("2021-12-31"), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red") +
  annotate("text", x = as.Date("2021-04-01"), y = Inf, label = "Periodo\nRetiros", 
           vjust = 1.5, hjust = 0.5, size = 3, color = "darkred", fontface = "bold") +
  scale_x_date(
    date_breaks = "2 years", 
    labels = function(x) ifelse(is.na(x), "", paste0(lubridate::year(x), "-T", lubridate::quarter(x)))
  )  +  scale_color_manual(values = c("W=60" = "gray50", "W=62" = "steelblue", "W=64" = "forestgreen")) +
  scale_linetype_manual(values = c("W=60" = "dashed", "W=62" = "solid", "W=64" = "dotdash")) +
  labs(title = "Total Connectivity Index (TCI)", x = "Periodo", y = "TCI (%)", color = "Ventana", linetype = "Ventana") +
  tema_graficos


# Grafico 2: NET
p_net <- ggplot(df_net_rw, aes(x = fecha, y = value, color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.6) + facet_wrap(~ sector, scales = "free_y", ncol = 2) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = as.Date("2008-09-01"), color = "gray40", linetype = "dotted") +
  geom_vline(xintercept = as.Date("2019-10-01"), color = "gray40", linetype = "dotted") +
  #retiros
  geom_vline(xintercept = as.Date("2020-07-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = as.Date("2020-12-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = as.Date("2021-04-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  annotate("rect", xmin = as.Date("2020-07-01"), xmax = as.Date("2021-12-31"), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red") +
  annotate("text", x = as.Date("2021-04-01"), y = Inf, label = "Periodo\nRetiros", 
           vjust = 1.5, hjust = 0.5, size = 2, color = "darkred", fontface = "bold") +
  scale_x_date(
    date_breaks = "2 years", 
    labels = function(x) ifelse(is.na(x), "", paste0(lubridate::year(x), "-T", lubridate::quarter(x)))
  )  +  scale_color_manual(values = c("W=60" = "gray50", "W=62" = "steelblue", "W=64" = "forestgreen")) +
  scale_linetype_manual(values = c("W=60" = "dashed", "W=62" = "solid", "W=64" = "dotdash")) +
  labs(title = "Net Directional Connectedness por Sector", subtitle = "Valores > 0: Transmisor Neto | Valores < 0: Receptor Neto",
       x = "Periodo", y = "Net Connectedness", color = "Ventana", linetype = "Ventana") +
  tema_graficos

# Grafico 3: TO
p_to <- ggplot(df_to_rw, aes(x = fecha, y = value, color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.6) + facet_wrap(~ sector, scales = "free_y", ncol = 2) +
  geom_vline(xintercept = as.Date("2008-09-01"), color = "gray40", linetype = "dotted") +
  geom_vline(xintercept = as.Date("2019-10-01"), color = "gray40", linetype = "dotted") +
  #retiros
  geom_vline(xintercept = as.Date("2020-07-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = as.Date("2020-12-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = as.Date("2021-04-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  annotate("rect", xmin = as.Date("2020-07-01"), xmax = as.Date("2021-12-31"), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red") +
  annotate("text", x = as.Date("2021-04-01"), y = Inf, label = "Periodo\nRetiros", 
           vjust = 1.5, hjust = 0.5, size = 3, color = "darkred", fontface = "bold") +
  scale_x_date(
    date_breaks = "2 years", 
    labels = function(x) ifelse(is.na(x), "", paste0(lubridate::year(x), "-T", lubridate::quarter(x)))
  )  +  scale_color_manual(values = c("W=60" = "gray50", "W=62" = "steelblue", "W=64" = "forestgreen")) +
  scale_linetype_manual(values = c("W=60" = "dashed", "W=62" = "solid", "W=64" = "dotdash")) +
  labs(title = "To-Connectedness por Sector (Intensidad Agregada de Transmision de Shocks)", subtitle = "Muestra qué proporcion de la varianza de otros sectores se explica por shocks de ese sector", 
       x = "Periodo", y = "To-Connectedness", color = "Ventana", linetype = "Ventana") +
  tema_graficos

# Grafico 4: FROM (EL NUEVO GRAFICO)
p_from <- ggplot(df_from_rw, aes(x = fecha, y = value, color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.6) + facet_wrap(~ sector, scales = "free_y", ncol = 2) +
  geom_vline(xintercept = as.Date("2008-09-01"), color = "gray40", linetype = "dotted") +
  geom_vline(xintercept = as.Date("2019-10-01"), color = "gray40", linetype = "dotted") +
  #retiros
  geom_vline(xintercept = as.Date("2020-07-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = as.Date("2020-12-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = as.Date("2021-04-01"), color = "red",
             linetype = "dotted", linewidth = 0.4) +
  annotate("rect", xmin = as.Date("2020-07-01"), xmax = as.Date("2021-12-31"), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red") +
  annotate("text", x = as.Date("2021-04-01"), y = Inf, label = "Periodo\nRetiros", 
           vjust = 1.5, hjust = 0.5, size = 3, color = "darkred", fontface = "bold") +
  scale_x_date(
    date_breaks = "2 years", 
    labels = function(x) ifelse(is.na(x), "", paste0(lubridate::year(x), "-T", lubridate::quarter(x)))
  )  +  scale_color_manual(values = c("W=60" = "gray50", "W=62" = "steelblue", "W=64" = "forestgreen")) +
  scale_linetype_manual(values = c("W=60" = "dashed", "W=62" = "solid", "W=64" = "dotdash")) +
  labs(
    title = "From-Connectedness por Sector (Recepcion de Shocks)",
    subtitle = "Muestra qué proporcion de la varianza del sector se explica por shocks de otros sectores.",
    x = "Periodo", y = "From-Connectedness (%)", color = "Ventana", linetype = "Ventana"
  ) +
  tema_graficos

# Guardar graficos
ggsave(file.path(ruta_output, "figs", "bcch_01_tci_rolling_w60.png"),  p_tci,  width = 10, height = 6, dpi = 150)
ggsave(file.path(ruta_output, "figs", "bcch_02_net_rolling_w60.png"),  p_net,  width = 10, height = 8, dpi = 150)
ggsave(file.path(ruta_output, "figs", "bcch_03_to_rolling_w60.png"),   p_to,   width = 10, height = 8, dpi = 150)
ggsave(file.path(ruta_output, "figs", "bcch_04_from_rolling_w60.png"), p_from, width = 10, height = 8, dpi = 150)

# Grafico 5: Pairwise Bancos vs FP
p_pairwise <- ggplot(df_pairwise_rw, aes(x = fecha, y = pairwise_value, color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ direccion, scales = "free_y", ncol = 1) + # Divide el grafico en 2 paneles (arriba y abajo)
  geom_vline(xintercept = as.Date("2008-09-01"), color = "gray40", linetype = "dotted") +
  geom_vline(xintercept = as.Date("2019-10-01"), color = "gray40", linetype = "dotted") +
  annotate("rect", xmin = as.Date("2020-07-01"), xmax = as.Date("2021-12-31"), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red") +
  annotate("text", x = as.Date("2021-04-01"), y = Inf, label = "Periodo\nRetiros", 
           vjust = 1.5, hjust = 0.5, size = 3, color = "darkred", fontface = "bold") +
  scale_x_date(
    date_breaks = "2 years", 
    labels = function(x) ifelse(is.na(x), "", paste0(lubridate::year(x), "-T", lubridate::quarter(x)))
  )  +  scale_color_manual(values = c("W=60" = "gray50", "W=62" = "steelblue", "W=64" = "forestgreen")) +
  scale_linetype_manual(values = c("W=60" = "dashed", "W=62" = "solid", "W=64" = "dotdash")) +
  labs(
    title = "Conectividad Bilateral (Pairwise): Bancos vs Fondos de Pensiones",
    subtitle = "Permite ver la direccion exacta del shock entre estos dos actores",
    x = "Periodo", y = "Varianza Explicada (%)", color = "Ventana", linetype = "Ventana"
  ) +
  tema_graficos # Reutiliza el tema que definimos en el paso anterior

ggsave(file.path(ruta_output, "figs", "bcch_05_pairwise_fp_bancos.png"), p_pairwise, width = 10, height = 8, dpi = 150)

# Grafico 6: Pairwise RestoMundo vs Bancos
p_pairwise_rm_bancos <- ggplot(df_pairwise_rm_Bancos, aes(x = fecha, y = pairwise_value, color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ direccion, scales = "free_y", ncol = 1) + # Divide el grafico en 2 paneles (arriba y abajo)
  geom_vline(xintercept = as.Date("2008-09-01"), color = "gray40", linetype = "dotted") +
  geom_vline(xintercept = as.Date("2019-10-01"), color = "gray40", linetype = "dotted") +
  annotate("rect", xmin = as.Date("2020-07-01"), xmax = as.Date("2021-12-31"), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red") +
  annotate("text", x = as.Date("2021-04-01"), y = Inf, label = "Periodo\nRetiros", 
           vjust = 1.5, hjust = 0.5, size = 3, color = "darkred", fontface = "bold") +
  scale_x_date(
    date_breaks = "2 years", 
    labels = function(x) ifelse(is.na(x), "", paste0(lubridate::year(x), "-T", lubridate::quarter(x)))
  )  +  scale_color_manual(values = c("W=60" = "gray50", "W=62" = "steelblue", "W=64" = "forestgreen")) +
  scale_linetype_manual(values = c("W=60" = "dashed", "W=62" = "solid", "W=64" = "dotdash")) +
  labs(
    title = "Conectividad Bilateral (Pairwise): Resto Mundo vs Bancos",
    subtitle = "Permite ver la direccion exacta del shock entre estos dos actores",
    x = "Periodo", y = "Varianza Explicada (%)", color = "Ventana", linetype = "Ventana"
  ) +
  tema_graficos # Reutiliza el tema que definimos en el paso anterior

ggsave(file.path(ruta_output, "figs", "bcch_06_pairwise_rm_bancos.png"), p_pairwise_rm_bancos, width = 10, height = 8, dpi = 150)

# Grafico 7: Pairwise Resto Mundo vs FP
p_pairwise_rm_FP <- ggplot(df_pairwise_rm_FP, aes(x = fecha, y = pairwise_value, color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ direccion, scales = "free_y", ncol = 1) + # Divide el grafico en 2 paneles (arriba y abajo)
  geom_vline(xintercept = as.Date("2008-09-01"), color = "gray40", linetype = "dotted") +
  geom_vline(xintercept = as.Date("2019-10-01"), color = "gray40", linetype = "dotted") +
  annotate("rect", xmin = as.Date("2020-07-01"), xmax = as.Date("2021-12-31"), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red") +
  annotate("text", x = as.Date("2021-04-01"), y = Inf, label = "Periodo\nRetiros", 
           vjust = 1.5, hjust = 0.5, size = 3, color = "darkred", fontface = "bold") +
  scale_x_date(
    date_breaks = "2 years", 
    labels = function(x) ifelse(is.na(x), "", paste0(lubridate::year(x), "-T", lubridate::quarter(x)))
  )  +  scale_color_manual(values = c("W=60" = "gray50", "W=62" = "steelblue", "W=64" = "forestgreen")) +
  scale_linetype_manual(values = c("W=60" = "dashed", "W=62" = "solid", "W=64" = "dotdash")) +
  labs(
    title = "Conectividad Bilateral (Pairwise): Resto Mundo vs Fondos de Pensiones",
    subtitle = "Permite ver la direccion exacta del shock entre estos dos actores",
    x = "Periodo", y = "Varianza Explicada (%)", color = "Ventana", linetype = "Ventana"
  ) +
  tema_graficos # Reutiliza el tema que definimos en el paso anterior

ggsave(filename = file.path(ruta_output, "figs", "bcch_07_pairwise_RM_FP.png"),      plot = p_pairwise_rm_FP,     width = 10, height = 8, dpi = 150, bg = "white")

# Grafico 8: Pairwise FP vs Banco Central
p_pairwise_fp_bcch <- ggplot(df_pairwise_fp_bcch, aes(x = fecha, y = pairwise_value, color = ventana, linetype = ventana)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ direccion, scales = "free_y", ncol = 1) +
  geom_vline(xintercept = as.Date("2008-09-01"), color = "gray40", linetype = "dotted") +
  geom_vline(xintercept = as.Date("2019-10-01"), color = "gray40", linetype = "dotted") +
  annotate("rect", xmin = as.Date("2020-07-01"), xmax = as.Date("2021-12-31"), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red") +
  annotate("text", x = as.Date("2021-04-01"), y = Inf, label = "Periodo\nRetiros", 
           vjust = 1.5, hjust = 0.5, size = 3, color = "darkred", fontface = "bold") +
  scale_x_date(
    date_breaks = "2 years", 
    labels = function(x) ifelse(is.na(x), "", paste0(lubridate::year(x), "-T", lubridate::quarter(x)))
  ) +
  scale_color_manual(values = c("W=60" = "gray50", "W=62" = "steelblue", "W=64" = "forestgreen")) +
  labs(
    title = "Conectividad Bilateral (Pairwise): Fondos de Pensiones vs Banco Central",
    subtitle = "Permite ver la direccion exacta del shock entre estos dos actores",
    x = "Periodo", y = "Varianza Explicada (%)", color = "Ventana", linetype = "Ventana"
  ) +
  tema_graficos

ggsave(file.path(ruta_output, "figs", "bcch_08_pairwise_fp_bcch.png"), 
       p_pairwise_fp_bcch, width = 10, height = 8, dpi = 150, bg = "white")

cat("\nPROCESO COMPLETADO EXITOSAMENTE. GRAFICOS 1 A 8 GENERADOS.\n")

# ============================================================
# SCRIPT R: Extraer valores puntuales en 4 fechas para tabla resumen
# Fuente: CSVs generados por 01_dy_pasos_1_a_5.R
# ============================================================
#install.packages("writexl")
library(tidyverse)
library(lubridate)
library(writexl)   # install.packages("writexl")

# Fechas de interés
fechas_target <- as.Date(c("2019-10-01","2020-04-01","2020-10-01","2022-01-01","2025-01-01"))
periodos_lbl  <- c("2019T4","2020T2","2020T4","2022T1","2025T1")

# 1. TCI (W=62)
tci <- read_csv(file.path(ruta_output,"bcch_05_tci_rolling_window.csv")) |>
  filter(ventana=="W=62", fecha %in% fechas_target) |>
  select(fecha, TCI_W62=tci)

# 2. From FP (W=62)
from_fp <- read_csv(file.path(ruta_output,"bcch_08_from_rolling_window.csv")) |>
  filter(ventana=="W=62", sector=="FP", fecha %in% fechas_target) |>
  select(fecha, From_FP_W62=value)

# 3. To FP (W=62)
to_fp <- read_csv(file.path(ruta_output,"bcch_07_to_rolling_window.csv")) |>
  filter(ventana=="W=62", sector=="FP", fecha %in% fechas_target) |>
  select(fecha, To_FP_W62=value)

# 4. Pairwise FP→Bancos
pw_fp_bancos <- read_csv(file.path(ruta_output,"bcch_09_pairwise_fp_bancos.csv")) |>
  filter(ventana=="W=62", direccion=="FP contagia a Bancos",
         fecha %in% fechas_target) |>
  select(fecha, Pairwise_FP_Bancos=pairwise_value)

# 5. Pairwise FP→RoW
pw_fp_row <- read_csv(file.path(ruta_output,"bcch_11_pairwise_RM_FP.csv")) |>
  filter(ventana=="W=62", direccion=="FP contagia a Resto Mundo",
         fecha %in% fechas_target) |>
  select(fecha, Pairwise_FP_RoW=pairwise_value)

# 6. Pairwise Bancos→FP
pw_bancos_fp <- read_csv(file.path(ruta_output,"bcch_09_pairwise_fp_bancos.csv")) |>
  filter(ventana=="W=62", direccion=="Bancos contagia a FP",
         fecha %in% fechas_target) |>
  select(fecha, Pairwise_Bancos_FP=pairwise_value)

# 7. Pairwise FP→BancoCentral
pw_fp_bcch <- read_csv(file.path(ruta_output, "bcch_12_pairwise_fp_bcch.csv")) |>
  filter(ventana == "W=62", direccion == "FP contagia a BancoCentral",
         fecha %in% fechas_target) |>
  select(fecha, Pairwise_FP_BCCh = pairwise_value)

# 8. Pairwise BancoCentral→FP
pw_bcch_fp <- read_csv(file.path(ruta_output, "bcch_12_pairwise_fp_bcch.csv")) |>
  filter(ventana == "W=62", direccion == "BancoCentral contagia a FP",
         fecha %in% fechas_target) |>
  select(fecha, Pairwise_BCCh_FP = pairwise_value)

# Unir todo
tabla <- list(tci, from_fp, to_fp, pw_fp_bancos, pw_fp_row, pw_bancos_fp,
              pw_fp_bcch, pw_bcch_fp) |>
  reduce(left_join, by = "fecha") |>
  mutate(
    periodo     = periodos_lbl[match(fecha, fechas_target)],
    Net_FP_W62  = To_FP_W62 - From_FP_W62,
    across(starts_with("Pairwise"), ~ .x)
  ) |>
  select(periodo, fecha,
         TCI_W62, From_FP_W62, To_FP_W62, Net_FP_W62,
         Pairwise_FP_Bancos, Pairwise_Bancos_FP,
         Pairwise_FP_RoW,
         Pairwise_FP_BCCh, Pairwise_BCCh_FP)

print(tabla)

# Guardar para pegar en Excel
write_csv(tabla, file.path(ruta_output, "bcch_W62_TABLA_RESUMEN_4FECHAS.csv"))
write_xlsx(tabla, file.path(ruta_output, "bcch_W62_TABLA_RESUMEN_4FECHAS.xlsx"))
cat("\n✓ Archivo guardado: TABLA_RESUMEN_4FECHAS.csv\n")
cat("\n--- VALORES PARA PEGAR EN EXCEL ---\n")