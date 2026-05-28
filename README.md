# Contagio sectorial y NBFIs en Chile

**Archivos de replicación** — Tesis para optar al grado de Magíster en Economía,
Pontificia Universidad Católica de Chile.

> **Título:** *Contagio sectorial y NBFIs en Chile: evidencia desde las Cuentas Nacionales y los retiros previsionales*
>
> **Autor:** Carlos Gonzalez
> **Profesores guías:** Javier Turen y Alejandro Jara
> **Año:** 2026
> **Instituto de Economía UC**

---

## Resumen
Este artículo analiza la estructura de intermediación financiera sectorial
en Chile, con énfasis en los intermediarios financieros no bancarios (NBFIs)
y en el contagio predictivo sectorial durante el episodio de retiros de fondos
previsionales: 2020T3-2021T4. Para ello, se recurre a las Cuentas Nacionales
por Sector Institucional del Banco Central de Chile, que reporta la matriz
who-to-whom de exposiciones bilaterales con frecuencia trimestral
para 2003T1-2025T3. 

Por medio del marco de contagio predictivo de Diebold y Yilmaz sobre la tasa
de crecimiento trimestral real del portafolio de intermediación de cada sector.
Se hallan tres regularidades. Primero, la conectividad agregada es estable en
promedio pero se duplica transitoriamente durante el período de retiros
previsionales y pandemia. Segundo, contrario a lo esperado, los Fondos de Pensiones
revierten su rol, pasando de transmisor neto en el período pre-retiros a receptor
neto, mientras el Banco Central emerge como principal transmisor sugiriendo un
rol relevante en sus medidas de mitigación capaz de revertir posiciones netas de contagio. 

---

## Estructura del repositorio

```
.
├── code/         # Scripts de Stata, R y Python para replicar resultados
|   ├── 2.1_do_files/         # Código en stata (grafico motivacion FSB)
│   ├── 2.2_ipynbs/           # Código en Python (limpieza de datos, gráficos y grafos)
│   └── 2.3_Rs/               # Código en R (modelo DY, estimaciones)
|
├── data/         # Datos
│   ├── 0_data_original/      # Datos crudos
│   └── 1_clean_data/         # Datos procesados, listos para análisis
|
├── output/       # Tablas y figuras generadas por el código
│   ├── 3_resultados_he_red/  # figuras, tablas y grafos
│   └── 3.2_resultado_DY/     # distintas salidas de metodología DY
|
├── 4_documentos/ # Documentos de tesis
│   ├── 4_tesis/              # Tesis formato pdf
│   └── 4.1_presentacion/     # Presentacion formato pdf
|
├── requirements.txt   # Dependencias de Python
├── renv.lock          # Dependencias de R (vía renv)
└── README.md
```

## Datos

Las fuentes utilizadas son:

- **Cuentas Nacionales por Sector Institucional (CNSI)** — Banco Central de Chile.
  Disponibles en: <https://www.bcentral.cl/>
- **Retiros previsionales** — [Fuente, p. ej. Superintendencia de Pensiones].

Para reproducir el analisis desde cero, descargar las series desde las fuentes
listadas y ubicarlas en `data/raw/` siguiendo la convención de nombres descrita
en `docs/data_sources.md`.

## Reproducción
...

### Requisitos

- **R** ≥ 4.3 (se usa `renv` para gestión de paquetes)
- **Python** ≥ 3.11

### Instalación

```bash
# Clonar el repo
git clone https://github.com/[tu-usuario]/replication-contagio-nbfi-chile.git
cd replication-contagio-nbfi-chile

# Entorno Python
python -m venv .venv
source .venv/bin/activate    # En Windows: .venv\Scripts\activate
pip install -r requirements.txt

# Entorno R — restaurar paquetes con renv
R -e "renv::restore()"
```

### Ejecución

Los scripts están numerados en orden de ejecución:

```bash
# 1. Limpieza y construcción de la base de datos (Python)
python code/python/Consolidado_data.ipynb

# 2. Análisis y estimaciones (Python y R)
  # 2.1 Generación de tablas y figuras (hechos estilizados y descriptivos)
    python code/python/MI_HE.ipynb
    python code/python/is_indice_red.ipynb
    python code/python/grafos_episodios_red.ipynb
    python code/python/plot_bonos.ipynb
# 2.2 Generación de tablas y figuras (método Diebold Yilmaz)
    Rscript code/R/02_bootstrap_dy.R #principal estimación DY
    Rscript code/R/01_dy_BCCh_distintasw # robustez con distintas ventanas
    Rscript code/R/01_dy_adaptive_lasso_bcch # robustez con lasso adaptative

Todos los resultados se generan en `output/` dentro de `3_resultados_he_red' o
`3.2_resultado_DY' según corresponda.

## Cita

Si utilizas estos materiales, por favor cita la tesis:

```bibtex
@mastersthesis{carlosgonzalez2026contagio,
  author  = {Carlos Gonzalez},
  title   = {Contagio sectorial y NBFIs en Chile: evidencia desde las Cuentas Nacionales y los retiros previsionales},
  school  = {Pontificia Universidad Católica de Chile, Instituto de Economía},
  year    = {2026},
  type    = {Tesis de Magíster en Economía}
}
```

## Licencia

El **código** se distribuye bajo licencia MIT (ver `LICENSE`).
Los **datos**, cuando aplique, están sujetos a los términos de las fuentes originales.

## Contacto

Carlos González - carlos.gonzlez@uc.cl
Profesores guías: 
      Javier Turén - Instituto de Economía UC
      Alejandro Jara - Banco Central de Chile
      

---

*Última actualización: 28 de mayo 2026.*
