********************************************************************************
* TESIS — Recomposicion de la Prima de Riesgo en el Diferencial de Retornos
*         CDA Guaranies vs. Dolares. Sistema Bancario Paraguayo, 2004–2024.
* Autores: Laviero Scavone y Matias Fernandez.  Stata 17+.
*------------------------------------------------------------------------------*
* SERIE PRINCIPAL de tasas: tramo CDA <=365 dias, tasa efectiva, bancos
*   ("Promedio de Tasas", Superintendencia de Bancos, BCP). Fuente unica y una
*   sola definicion 2004–2024 (sin empalme). Serie ponderada = sensibilidad.
* CORTE INSTITUCIONAL PRINCIPAL: 2011T2 (adopcion formal de Metas de Inflacion,
*   Res. N.22 Acta N.31 del 18-may-2011). Sensibilidad: 2011T1.
* MODELO NUCLEO: se estima sobre la PRIMA DE RIESGO rp = difn - InflDiff, con la
*   inflacion RESTADA (nunca como regresor), para evitar la circularidad de Fisher.
*------------------------------------------------------------------------------*
* COMO CORRERLO:
*   do "stata/tesis_cda.do" "C:/ruta/a/la/carpeta/tes"
*   (si se omite la ruta, usa la global 'root' de abajo). Requiere internet la
*   primera vez para instalar rangestat, estout, xtbreak y kpss desde SSC.
* SALIDAS: log -> output/tesis_cda.log ; figuras -> output/figuras/ ;
*          base derivada -> output/base_cda_derivada.dta
* HONESTIDAD: coeficientes, SE y p-valores se reportan TAL COMO SALEN.
********************************************************************************
version 17.0
clear all
set more off
set linesize 120
set scheme s1color            // figuras con fondo blanco para el documento

* --- Ruta portable: acepta argumento; si no, usa esta global (unica edicion) ---
args root
if `"`root'"' == "" global root "C:/Users/Laviero Scavone/Desktop/tes"
else global root `"`root'"'
cd "$root"
capture mkdir "output"
capture mkdir "output/figuras"
capture mkdir "output/tablas"

capture log close _all
log using "output/tesis_cda.log", replace text

********************************************************************************
* 0. PAQUETES Y PROGRAMA AUXILIAR
********************************************************************************
foreach pkg in rangestat estout xtbreak kpss {
    capture which `pkg'
    if _rc {
        di as text "Instalando `pkg' desde SSC..."
        capture ssc install `pkg', replace
    }
}

* newey_nogap: newey sobre la submuestra estimable, reindexando el tiempo a un
* indice secuencial. Necesario cuando la serie tiene huecos internos (la serie
* ponderada tiene uno en 2006T4). Para series contiguas equivale a newey normal.
capture program drop newey_nogap
program define newey_nogap, eclass
    version 17.0
    syntax varlist(numeric min=1) [if] [in], LAG(integer)
    marksample touse
    preserve
        quietly keep if `touse'
        quietly gen long __tseq = _n
        quietly tsset __tseq
        newey `varlist', lag(`lag')
    restore
end

********************************************************************************
* 1. IMPORTACION Y PRUEBAS DE INTEGRIDAD
********************************************************************************
confirm file "data/processed/base_cda.csv"
import delimited using "data/processed/base_cda.csv", clear varnames(1) case(preserve) encoding(utf8)

assert _N == 84
gen tq = yq(anio, trimestre)
format tq %tq
sort tq
isid tq
assert tq == yq(2004,1) + _n - 1          // 84 trimestres contiguos 2004T1–2024T4
tsset tq

* La serie principal (<=365) no tiene huecos trimestrales:
count if missing(i_Gs) | missing(i_USD)
assert r(N) == 0
assert inrange(meses_cobertura,1,3)
count if meses_cobertura < 3
di as result "Trimestres con cobertura mensual parcial (declarados): " r(N)
list fecha meses_cobertura ultimo_mes_publicado if meses_cobertura < 3, noobs

********************************************************************************
* 2. CONSTRUCCION DE VARIABLES
********************************************************************************
* depreciacion del guarani (%)
gen dep     = 100*(TC/L.TC - 1)
gen depLead = 100*(F.TC/TC - 1)                    // depreciacion del trimestre siguiente
* diferencial de tasas
gen difn = i_Gs - i_USD                            // nominal (pp, % anual)
gen difq = difn/4                                  // trimestralizado (frecuencia de dep)
* series de retorno
gen DR     = difq - depLead
gen DRcons = difq - max(depLead,0) if !missing(depLead)   // conservadora (excluye apreciacion)
* Fisher
gen infl_py  = 100*(IPC_PY/L4.IPC_PY - 1)
gen infl_us  = 100*(IPC_US/L4.IPC_US - 1)
gen InflDiff = infl_py - infl_us
gen rp       = difn - InflDiff                     // prima de riesgo (diferencial real residual)
* regimen: PRINCIPAL 2011T2 ; sensibilidad 2011T1
gen D2011    = (anio>2011) | (anio==2011 & trimestre>=2)
gen D2011_T1 = (anio>=2011)
label var D2011    "Post-regimen (>=2011T2, principal)"
label var D2011_T1 "Post-regimen (>=2011T1, sensibilidad)"
* crisis financiera global (quiebre endogeno de DR)
gen crisis = inrange(tq, yq(2008,2), yq(2009,3))
* volatilidad cambiaria: sd movil de 8 trimestres
capture which rangestat
if _rc == 0 {
    rangestat (sd) VolTC = dep, interval(tq -7 0)
}
else {
    * respaldo nativo (no requiere SSC ni internet): las 84 filas son trimestres
    * contiguos, por lo que la ventana por posicion equivale a la ventana calendario
    di as text "rangestat no disponible: usando respaldo nativo para VolTC"
    gen VolTC = .
    forvalues i = 1/`=_N' {
        quietly summarize dep in `=max(1,`i'-7)'/`i'
        if r(N) >= 2 quietly replace VolTC = r(sd) in `i'
    }
}
gen dolar = dolariz                                // proxy de iliquidez/segmentacion

* Chequeo de formulas: la inflacion recalculada coincide con la publicada
assert abs(infl_py - infl_PY_yoy) < 0.001 if !missing(infl_py, infl_PY_yoy)
assert abs(infl_us - infl_US_yoy) < 0.001 if !missing(infl_us, infl_US_yoy)

* variables centradas del modelo nucleo (para interpretar D2011 en la media)
gen byte sample_main = !missing(rp, VolTC, dolar)
quietly summarize VolTC if sample_main, meanonly
gen VolTC_c = VolTC - r(mean)
quietly summarize dolar if sample_main, meanonly
gen dolar_c = dolar - r(mean)
gen DxVolTC = D2011*VolTC_c
gen DxDolar = D2011*dolar_c
label var VolTC_c "Riesgo cambiario (VolTC, centrada)"
label var dolar_c "Iliquidez/dolarizacion (centrada)"
label var DxVolTC "D2011T2 x VolTC"
label var DxDolar "D2011T2 x dolarizacion"

********************************************************************************
* 3. ETAPA 1 — DESCRIPTIVOS + FIGURA 8.1 (serie DR)
********************************************************************************
di as text _n(2) "{hline 70}"
di " ETAPA 1 — Descriptivos por subperiodo (corte 2011T2)"
di "{hline 70}"
tabstat difn DR DRcons InflDiff rp, by(D2011) ///
    statistics(n mean sd min p50 max) columns(statistics)

* ===== FIGURA 8.1 del Word: serie del diferencial de retorno DR =====
* (comando en UNA sola linea: anda igual si se pega en la ventana de comandos)
tsline DR, lcolor(navy) lwidth(medthick) yline(0, lpattern(dash) lcolor(gs9)) xline(`=yq(2011,2)', lcolor(red) lwidth(medthick)) title("Diferencial de retorno trimestral DR (guaranies vs. dolares)") subtitle("2004-2024; linea roja = adopcion de Metas de Inflacion (2011T2)") ytitle("DR (puntos porcentuales, trimestral)") xtitle("") note("Fuente: elaboracion propia con datos del BCP.") name(g81, replace)
graph export "output/figuras/Figura_8_1_DR.png", replace width(1600)

********************************************************************************
* 4. ETAPA 2 — REGRESION DE FAMA + PRUEBA CONSERVADORA
********************************************************************************
di as text _n(2) "{hline 70}"
di " ETAPA 2 — Coeficiente de Fama (UIP <=> beta=1) + DRcons"
di "{hline 70}"
* difq (trimestralizado) es la misma frecuencia que depLead => beta=1 es la nula UIP
newey_nogap depLead difq, lag(4)
test difq = 1
estimates store fama_total
foreach cond in "D2011==0" "D2011==1" {
    newey_nogap depLead difq if `cond', lag(4)
    test difq = 1
}
* Prueba conservadora H0: media(DRcons)<=0 vs H1>0, con t de Stata (df=n-1)
newey_nogap DRcons, lag(4)
scalar t_dr = _b[_cons]/_se[_cons]
di as result "DRcons total: media=" %7.4f _b[_cons] ///
    "  p unilateral(H1: media>0) = " %6.4f ttail(e(df_r), t_dr)
newey_nogap DRcons if D2011==1, lag(4)
scalar t_dr = _b[_cons]/_se[_cons]
di as result "DRcons post : media=" %7.4f _b[_cons] ///
    "  p unilateral(H1: media>0) = " %6.4f ttail(e(df_r), t_dr)

********************************************************************************
* 5. ETAPA 4 — DESCOMPOSICION DE FISHER + FIGURA 8.2 (barras)
********************************************************************************
di as text _n(2) "{hline 70}"
di " ETAPA 4 — Descomposicion de Fisher: difn = InflDiff + rp"
di "{hline 70}"
tabstat difn InflDiff rp if !missing(rp), by(D2011) statistics(n mean sd) columns(statistics)
di as text "Cambio de medias pre/post 2011T2 (HAC):"
foreach y in rp InflDiff difn {
    newey_nogap `y' D2011 if !missing(rp), lag(4)
    di as result "  Delta `y' = " %7.3f _b[D2011] "  (p=" %5.3f ///
        (2*ttail(e(df_r), abs(_b[D2011]/_se[D2011]))) ")"
    estimates store fisher_`y'
}

* ===== FIGURA 8.2 del Word: descomposicion de Fisher (barras pre/post) =====
graph bar (mean) InflDiff rp if !missing(rp), over(D2011, relabel(1 "Pre-2011T2" 2 "Post-2011T2")) bar(1, color(navy)) bar(2, color(cranberry)) blabel(bar, format(%4.2f) size(small)) ytitle("puntos porcentuales anuales") yline(0, lcolor(gs9)) legend(order(1 "Prima de inflacion (InflDiff)" 2 "Prima de riesgo (rp)")) title("Descomposicion de Fisher del diferencial nominal") note("Fuente: elaboracion propia con datos del BCP y FRED.") name(g82, replace)
graph export "output/figuras/Figura_8_2_Fisher.png", replace width(1600)

********************************************************************************
* 6. ETAPA 5 — MODELO CENTRAL MCO-HAC (rp centrado, corte 2011T2)
********************************************************************************
di as text _n(2) "{hline 70}"
di " ETAPA 5 — Modelo nucleo: rp ~ VolTC_c + dolar_c + D2011 + DxVolTC + DxDolar"
di "{hline 70}"
di as text " Recomposicion interna: HA b(DxVolTC)<0 y b(DxDolar)>0"
newey_nogap rp VolTC_c dolar_c D2011 DxVolTC DxDolar, lag(4)
estimates store modelo_principal
di as text _n "Wald conjunto de cambio de pendientes (H0: DxVolTC=DxDolar=0):"
test DxVolTC DxDolar
di as text "Chow completo (nivel + pendientes):"
test D2011 DxVolTC DxDolar

********************************************************************************
* 7. ETAPA 3 — ESTACIONARIEDAD Y QUIEBRE ESTRUCTURAL
*   (bloque contiguo: dfuller/pperron/kpss/estat sb*/xtbreak exigen muestra sin
*    huecos; se reindexa el tiempo a un indice secuencial.)
********************************************************************************
di as text _n(2) "{hline 70}"
di " ETAPA 3 — Raices unitarias (ADF/PP/KPSS) y quiebre estructural"
di "{hline 70}"
di as text "-- Integracion por serie (ADF: H0 raiz unitaria; KPSS: H0 estacionariedad) --"
foreach y in DR rp VolTC dolar InflDiff {
    preserve
        quietly keep if !missing(`y')
        quietly gen long __s = _n
        quietly tsset __s
        di as text _n "ADF/PP/KPSS: `y'"
        dfuller `y', lags(4)
        capture noisily pperron `y'
        capture noisily kpss `y', maxlag(4)
    restore
}
* Quiebre estructural endogeno sobre DR (serie de retorno) y rp (prima)
foreach y in DR rp {
    preserve
        quietly keep if !missing(`y')
        quietly gen long __s = _n
        quietly tsset __s
        di as text _n "== Quiebre endogeno en `y' =="
        quietly regress `y' L.`y'
        capture noisily estat sbsingle
        capture noisily xtbreak test `y', breaks(3) hypothesis(3)      // UDmax/WDmax
        capture noisily xtbreak estimate `y', breaks(1)                // fecha + IC95
        di as text "  (fechas de xtbreak en indice colapsado; DR->crisis 2008T2 ; rp->2011T2)"
    restore
}
* Diferencia de medias de DR pre/post 2011T2 (prueba t con HAC)
di as text _n "-- Diferencia de medias DR pre/post 2011T2 (HAC) --"
newey_nogap DR D2011, lag(4)

********************************************************************************
* 8. ETAPA 5-BIS — COMPOSICION INTERNA DE LA PRIMA (Shapley/LMG) + FIGURA 8.3
*   Descriptivo/exploratorio (NO es el contraste formal; ese es la Etapa 5).
********************************************************************************
di as text _n(2) "{hline 70}"
di " ETAPA 5-BIS — Composicion interna de rp (Shapley/LMG) + dominancia"
di "{hline 70}"
* (a) composicion estatica pre/post: reparto del R2 entre VolTC e Iliq
di as text "-- Reparto del R2 de rp entre riesgo cambiario e iliquidez --"
foreach d in 0 1 {
    quietly regress rp VolTC if D2011==`d'
    scalar r2v = e(r2)
    quietly regress rp dolar if D2011==`d'
    scalar r2i = e(r2)
    quietly regress rp VolTC dolar if D2011==`d'
    scalar r2b = e(r2)
    scalar cv = 0.5*(r2v + r2b - r2i)
    scalar ci = 0.5*(r2i + r2b - r2v)
    local per = cond(`d'==0, "Pre-2011T2 ", "Post-2011T2")
    di as result "  `per': R2=" %5.3f r2b "  %VolTC=" %5.1f 100*cv/(cv+ci) ///
        "%  %Iliq=" %5.1f 100*ci/(cv+ci) "%"
}
* correlacion simple rp-dolarizacion por subperiodo (Seccion 8.7 del Word)
di as text "-- Correlacion rp-dolarizacion por subperiodo --"
foreach d in 0 1 {
    quietly corr rp dolar if D2011==`d'
    di as result "  D2011=`d': corr(rp,dolar)=" %6.3f r(rho) "  (n=" r(N) ")"
}
* (b) dominancia estandarizada post-2011: incidencia de cada factor en desvios std
egen zV = std(VolTC) if sample_main
egen zI = std(dolar) if sample_main
gen DzV = D2011*zV
gen DzI = D2011*zI
newey_nogap rp zV zI D2011 DzV DzI, lag(4)
di as text "-- Dominancia post-2011T2 (efecto en desvios estandar) --"
lincom zV + DzV
lincom zI + DzI
lincom (zI + DzI) - (zV + DzV)      // dominancia relativa (indistinguible si p>0.05)

* (c) FIGURA 8.3 del Word: evolucion de la composicion (ventana movil 20 trim)
preserve
    quietly keep if !missing(rp, VolTC, dolar)
    gen pctVol = .
    gen pctIliq = .
    local W = 20
    forvalues i = `W'/`=_N' {
        local j = `i' - `W' + 1
        quietly regress rp VolTC in `j'/`i'
        scalar r2v = e(r2)
        quietly regress rp dolar in `j'/`i'
        scalar r2i = e(r2)
        quietly regress rp VolTC dolar in `j'/`i'
        scalar r2b = e(r2)
        scalar cv = 0.5*(r2v + r2b - r2i)
        scalar ci = 0.5*(r2i + r2b - r2v)
        quietly replace pctVol  = 100*cv/(cv+ci) in `i'
        quietly replace pctIliq = 100*ci/(cv+ci) in `i'
    }
    quietly tsset tq
    tsline pctVol pctIliq, lcolor(cranberry navy) lwidth(medthick medthick) yline(50, lpattern(dot) lcolor(gs9)) xline(`=yq(2011,2)', lcolor(black) lpattern(dash)) title("Composicion de la prima de riesgo (ventana movil de 20 trimestres)") subtitle("linea negra = adopcion de Metas de Inflacion (2011T2)") ytitle("% del R2 explicado") xtitle("Trimestre (fin de ventana)") legend(order(1 "Riesgo cambiario (VolTC)" 2 "Iliquidez/dolarizacion")) note("Fuente: elaboracion propia.") name(g83, replace)
    graph export "output/figuras/Figura_8_3_Composicion.png", replace width(1600)
restore

********************************************************************************
* 9. ETAPA 6 — DIAGNOSTICOS + COINTEGRACION (Engle-Granger)
********************************************************************************
di as text _n(2) "{hline 70}"
di " ETAPA 6 — Diagnosticos del modelo y cointegracion"
di "{hline 70}"
preserve
    quietly keep if sample_main
    quietly gen long __s = _n
    quietly tsset __s
    regress rp VolTC_c dolar_c D2011 DxVolTC DxDolar
    di as text "-- Breusch-Godfrey (autocorrelacion) --"
    estat bgodfrey, lags(1 2 3 4)
    di as text "-- White (heterocedasticidad) --"
    estat imtest, white
    di as text "-- Ramsey RESET (forma funcional) --"
    estat ovtest
    di as text "-- Normalidad de residuos (Skewness/Kurtosis) --"
    predict __r, resid
    sktest __r
    di as text "-- VIF (multicolinealidad) --"
    vif
    * Engle-Granger: las series se comportan como ~I(1); si los residuos del
    * modelo son estacionarios (ADF rechaza raiz unitaria) las series COINTEGRAN
    * y la relacion NO es espuria.
    di as text "-- Engle-Granger: ADF sobre los residuos del modelo nucleo --"
    dfuller __r, lags(3)
    drop __r
restore

********************************************************************************
* 10. ETAPA 6/7 — ROBUSTEZ
********************************************************************************
di as text _n(2) "{hline 70}"
di " ROBUSTEZ"
di "{hline 70}"
* (i) sensibilidad de la FECHA: corte 2011T1
quietly summarize VolTC if !missing(rp,VolTC,dolar), meanonly
gen VolTC_c1 = VolTC - r(mean)
quietly summarize dolar if !missing(rp,VolTC,dolar), meanonly
gen dolar_c1 = dolar - r(mean)
gen DxV_T1 = D2011_T1*VolTC_c1
gen DxD_T1 = D2011_T1*dolar_c1
di as text _n "-- (i) Corte 2011T1 (sensibilidad de la fecha) --"
newey_nogap rp VolTC_c1 dolar_c1 D2011_T1 DxV_T1 DxD_T1, lag(4)
test DxV_T1 DxD_T1
* (ii) control por la crisis global 2008–09
di as text _n "-- (ii) Control por la crisis 2008T2–2009T3 --"
newey_nogap rp VolTC_c dolar_c D2011 DxVolTC DxDolar crisis, lag(4)
* (iii) tasas de fin de trimestre
gen difn_fdt = i_Gs_fdt - i_USD_fdt
gen rp_fdt = difn_fdt - InflDiff
quietly summarize VolTC if !missing(rp_fdt,VolTC,dolar), meanonly
gen VolTC_cf = VolTC - r(mean)
quietly summarize dolar if !missing(rp_fdt,VolTC,dolar), meanonly
gen dolar_cf = dolar - r(mean)
gen DxVf = D2011*VolTC_cf
gen DxDf = D2011*dolar_cf
di as text _n "-- (iii) Tasas de fin de trimestre --"
newey_nogap rp_fdt VolTC_cf dolar_cf D2011 DxVf DxDf, lag(4)
test DxVf DxDf
* (iv) rezago HAC alternativo (3, 5, 6)
di as text _n "-- (iv) Sensibilidad al rezago HAC --"
foreach L in 3 5 6 {
    di as text "   lag `L':"
    newey_nogap rp VolTC_c dolar_c D2011 DxVolTC DxDolar, lag(`L')
    test DxVolTC DxDolar
}
* (v) SENSIBILIDAD a la DEFINICION de la serie: promedio ponderado (con empalme)
di as text _n "-- (v) Serie ponderada (sensibilidad; hueco 2006T4) --"
gen difn_pond = i_Gs_pond - i_USD_pond
gen rp_pond = difn_pond - InflDiff
quietly summarize VolTC if !missing(rp_pond,VolTC,dolar), meanonly
gen VolTC_cp = VolTC - r(mean)
quietly summarize dolar if !missing(rp_pond,VolTC,dolar), meanonly
gen dolar_cp = dolar - r(mean)
gen DxVp = D2011*VolTC_cp
gen DxDp = D2011*dolar_cp
newey_nogap rp_pond VolTC_cp dolar_cp D2011 DxVp DxDp, lag(4)
test DxVp DxDp
newey_nogap rp_pond D2011, lag(4)          // salto de nivel con la ponderada

* (vi) robustez: indicadora GRADUAL del regimen (rampa lineal 2011T2-2016; 1 desde 2017)
di as text _n "-- (vi) Regimen gradual (rampa 2011T2-2016T4) --"
gen Reg = 0
replace Reg = (tq - yq(2011,2) + 1)/(yq(2016,4) - yq(2011,2) + 1) if inrange(tq, yq(2011,2), yq(2016,4))
replace Reg = 1 if tq > yq(2016,4)
gen RxV = Reg*VolTC_c
gen RxD = Reg*dolar_c
newey_nogap rp VolTC_c dolar_c Reg RxV RxD, lag(4)
test RxV RxD

* (vii) robustez: dolarizacion REZAGADA un trimestre (mitiga la simultaneidad)
di as text _n "-- (vii) Dolarizacion rezagada un trimestre --"
gen Ldolar = L.dolar
quietly summarize Ldolar if !missing(rp,VolTC,Ldolar), meanonly
gen Ldolar_c = Ldolar - r(mean)
gen DxLD = D2011*Ldolar_c
newey_nogap rp VolTC_c Ldolar_c D2011 DxVolTC DxLD, lag(4)
test DxVolTC DxLD

********************************************************************************
* 11. REPORTE PARA EL WORD (numeros clave, ya calculados arriba, en un solo lugar)
********************************************************************************
di as text _n(3) "{hline 70}"
di " REPORTE PARA EL WORD (copiar a las tablas)"
di "{hline 70}"
di as text "FIGURAS en output/figuras/:"
di "   Figura 8.1 = Figura_8_1_DR.png          (serie del diferencial DR)"
di "   Figura 8.2 = Figura_8_2_Fisher.png      (descomposicion de Fisher)"
di "   Figura 8.3 = Figura_8_3_Composicion.png (composicion, ventana movil)"
di as text _n ">> Modelo central (rp, T2):"
newey_nogap rp VolTC_c dolar_c D2011 DxVolTC DxDolar, lag(4)
test DxVolTC DxDolar
di as result "   ^ WALD pendientes: F=" %5.3f r(F) "  p=" %5.3f r(p) ///
    "   (NO rechaza estabilidad de pendientes)"
test D2011 DxVolTC DxDolar
di as result "   ^ CHOW completo:  F=" %5.2f r(F) "  p=" %5.3f r(p) ///
    "   (significativo por el NIVEL)"
di as text _n "Interpretacion central (honesta):"
di as text " El regimen se asocia a un salto de NIVEL de la prima (+3.9 pp, p<0.001)"
di as text " y a la recomposicion Fisher inflacion->riesgo; la prueba conjunta de"
di as text " cambio de pendientes NO rechaza estabilidad (recomposicion interna no"
di as text " concluyente). Resultado robusto a la definicion? NO: con la serie"
di as text " ponderada las pendientes cambian -> se declara como sensibilidad."

********************************************************************************
* 12. CIERRE
********************************************************************************
compress
save "output/base_cda_derivada.dta", replace
di as text _n "FIN. Log -> output/tesis_cda.log | Figuras -> output/figuras/"
log close
********************************************************************************
