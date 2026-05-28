*******************************************************
* algoritmo_nbfi_jurisdiccion_fsb.do
* Construye:
* (i) gráfico de barras agrupadas
* Datos: obtenidos del Global Monitoring Report on
* Nonbank Financial Intermediation: Data dashboard del FSB.
*******************************************************

clear all
set more off

*******************************************************
* 1) RUTA DE SALIDA
*******************************************************
local outgraph "../../3_resultados/3_resultados_he_red/figuras/1_MI_HE"


********************************************************************************
* Gráfico de barras agrupadas:
* Participación de NBFIs en activos financieros totales por jurisdicción (2024)
********************************************************************************

*--------------------------------------------------
* 2. Datos
*--------------------------------------------------
input str4 categoria Chile ADV EME
"NBFI" 54.6 57.9 27.5
"PF"   23.7 10.9  1.6
"IC"    8.4  8.4  5.9
"OFIs" 22.4 37.8 19.8
end

*--------------------------------------------------
* 3. Variable de orden para forzar:
*    NBFI, PF, IC, OFIs
*--------------------------------------------------
gen orden = .
replace orden = 1 if categoria == "NBFI"
replace orden = 2 if categoria == "PF"
replace orden = 3 if categoria == "IC"
replace orden = 4 if categoria == "OFIs"

label define orden_lbl 1 "NBFI" 2 "PF" 3 "IC" 4 "OFIs"
label values orden orden_lbl

sort orden

*--------------------------------------------------
* 4. Estilo general
*--------------------------------------------------
graph set window fontface "Times New Roman"

*--------------------------------------------------
* 5. Gráfico
*--------------------------------------------------
graph bar (asis) Chile ADV EME, ///
    over(orden, relabel(1 "NBFI" 2 "PF" 3 "IC" 4 "OFIs") label(labsize(large))) ///
    asyvars ///
    bargap(10) ///
    bar(1, color("46 84 134")   lcolor("46 84 134")) ///
    bar(2, color("230 126 34")  lcolor("230 126 34")) ///
    bar(3, color("160 160 160") lcolor("160 160 160")) ///
    ylabel(0(10)70, ///
        labsize(small) ///
        angle(0) ///
        grid ///
        glcolor(gs12) ///
        glwidth(medthin)) ///
    yscale(range(0 72)) ///
    ytitle("Porcentaje (%)", size(vsmed)) ///
     title("Participación en el total de activos financieros", size(medium) color(black)) ///
    legend(order(1 "Chile" 2 "ADV" 3 "EME") ///
        rows(1) ///
        position(6) ///
        ring(1) ///
        size(medium) ///
        region(lcolor(none) fcolor(none))) ///
    graphregion(color(white)) ///
    plotregion(color(white)) ///
    scheme(s1color) ///
	note("Fuente: elaboración propia en base a datos de FSB (2025).", size(vsmall))

graph export "`outgraph'/fig1_activos_jurisdicciones.png", replace width(3000)
