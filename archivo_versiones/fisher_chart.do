*==============================================================================*
* Etapa 4 - Descomposicion de Fisher del diferencial nominal (barras agrupadas)
* Correr con:   do "C:/Users/Laviero Scavone/Desktop/tes/stata/fisher_chart.do"
* Nota: se usan barras AGRUPADAS (no apiladas con %) porque la prima de riesgo
* pre-2011 es NEGATIVA (-1,41), y un apilado con porcentajes no tiene sentido.
*==============================================================================*
import delimited "C:/Users/Laviero Scavone/Desktop/tes/data/processed/base_cda.csv", clear case(preserve)

* --- variables derivadas (como en tesis_cda.do) ---
gen difn     = i_Gs - i_USD
gen InflDiff = infl_PY_yoy - infl_US_yoy
gen rp       = difn - InflDiff
gen D2011    = (anio >= 2011)
label define reg 0 "Pre-2011" 1 "Post-2011", replace
label values D2011 reg

* --- carpeta de salida (por si no existe) ---
cap mkdir "C:/Users/Laviero Scavone/Desktop/tes/output"
cap mkdir "C:/Users/Laviero Scavone/Desktop/tes/output/figuras"

graph bar (mean) InflDiff rp, over(D2011) ///
    bar(1, color(navy)) bar(2, color(orange)) ///
    blabel(bar, format(%4.2f)) ///
    yline(0, lcolor(gs8)) ///
    ytitle("Puntos porcentuales (% anual)") ///
    title("Etapa 4 - Descomposicion de Fisher del diferencial nominal") ///
    subtitle("Prima de inflacion vs prima de riesgo (medias por subperiodo)") ///
    legend(order(1 "Prima de inflacion (InflDiff)" 2 "Prima de riesgo (rp)") rows(1)) ///
    graphregion(color(white))

graph export "C:/Users/Laviero Scavone/Desktop/tes/output/figuras/etapa4_fisher_stata.png", replace width(1400)
