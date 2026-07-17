*==============================================================================*
*  TESIS: Cambio de Nivel y Composicion del Diferencial de Tasas entre CDA en   *
*         Guaranies y Dolares durante la Transicion al Regimen de Metas de      *
*         Inflacion: Paraguay, 2004-2024.   (Scavone & Fernandez)               *
*==============================================================================*
*  DO-FILE CANONICO - STATA 17.   Version: v2026-07-16  (post-auditoria integral)
*
*  UNA SOLA BASE, UNA SOLA CORRIDA, UNA SOLA SALIDA:
*   - Base canonica  : data/processed/base_cda.csv  (84 trimestres, 2004Q1-2024Q4)
*   - Serie primaria : tramo CDA <=365 dias (i_Gs, i_USD) - fuente unica SB/BCP,
*                      sin empalme. La ponderada (i_*_pond) es SOLO sensibilidad.
*   - Nucleo         : rp = (i_Gs - i_USD) - InflDiff  ("diferencial real
*                      residual", proxy de la prima de riesgo; NO prima pura:
*                      mezcla error de expectativas y primas de liquidez/credito).
*   - Regimen        : D2011 = 1 desde 2011:T2 (Resolucion N.22, Acta N.31 del
*                      18-may-2011). PRINCIPAL = T2; sensibilidad = T1.
*   - Cada bloque lleva el numero de la tabla del documento (Tabla 8.x) para
*     correspondencia 1:1 entre Word, do-file y salida.
*
*  VEREDICTO QUE REPRODUCE ESTE ARCHIVO (unico, sin ambiguedad):
*   (i)  El nivel del diferencial real residual SUBE +3,89 pp tras 2011T2
*        (p<0,001), robusto a todas las sensibilidades.
*   (ii) Las PENDIENTES de VolTC e Iliq NO cambian conjuntamente al 5%
*        (Wald F=2,00; p=0,142). NO hay evidencia de recomposicion interna.
*   (iii) El quiebre endogeno (Bai-Perron) de rp cae en 2011Q2.
*   Con la serie ponderada (sensibilidad) las pendientes SI cambian (p<0,001):
*   la evidencia de pendientes es sensible a la definicion de la serie; la de
*   nivel es robusta. Ninguna conclusion del documento excede estos limites.
*
*  COMO CORRERLO (portatil, sin rutas personales):
*    1) Abrir Stata 17 y situarse en la CARPETA RAIZ del proyecto:
*         cd "<carpeta-del-proyecto>"        (la carpeta que contiene data/ y stata/)
*    2) do "stata/tesis_cda.do"
*    Los modulos SSC (rangestat, estout, xtbreak, kpss) se instalan solos si faltan.
*    Toda la sesion queda registrada en output/log_tesis_cda.log
*==============================================================================*

clear all
set more off
version 17

*------------------------------------------------------------------------------*
* 0. VALIDACION DE CARPETA (portatil: NO hay rutas absolutas en este archivo)
*------------------------------------------------------------------------------*
capture confirm file "data/processed/base_cda.csv"
if _rc {
    di as error "ERROR: no se encuentra data/processed/base_cda.csv"
    di as error "Ejecute primero:  cd <carpeta raiz del proyecto>  (la que contiene data/ y stata/)"
    exit 601
}
capture mkdir "output"
capture mkdir "output/tablas"
capture mkdir "output/figuras"

* --- LOG de toda la sesion (evidencia de la corrida) ---
capture log close
log using "output/log_tesis_cda.log", replace text

di "Stata: `c(stata_version)'  |  Fecha de corrida: `c(current_date)' `c(current_time)'"

* --- Modulos SSC: instalar solo si faltan ---
foreach pkg in rangestat estout xtbreak kpss {
    capture which `pkg'
    if _rc capture ssc install `pkg', replace
    capture which `pkg'
    if _rc di as text "[aviso] `pkg' no disponible (sin internet?): los bloques que lo usan se omiten con aviso."
}

*------------------------------------------------------------------------------*
* Bandwidth HAC (Newey-West): lag = 4 (~un anio de memoria trimestral).
* Andrews O(T^{1/3}) ~ 4,4 con T=84; regla NW floor(4*(84/100)^{2/9}) = 3.
* Sensibilidad con lags 3/5/6 en el bloque 8.8.
*------------------------------------------------------------------------------*
global HAC = 4

*------------------------------------------------------------------------------*
* newey_nogap: newey sobre la submuestra estimable colapsada a indice secuencial.
*   newey exige muestra sin huecos temporales (r(498)). La SERIE PRIMARIA no
*   tiene huecos internos (la muestra del modelo es contigua): el colapso solo
*   remueve bordes (leads/lags) y NO altera la adyacencia temporal. En la serie
*   PONDERADA (sensibilidad, hueco 2006Q4) el colapso vuelve adyacentes 2006Q3
*   y 2007Q1: se declara como limitacion de ESA sensibilidad.
*   Es el mismo tratamiento que statsmodels aplica al vector contiguo sin NaN.
*------------------------------------------------------------------------------*
capture program drop newey_nogap
program define newey_nogap, eclass
    syntax varlist(numeric min=1) [if] [in], lag(integer 4)
    marksample touse
    preserve
    quietly keep if `touse'
    quietly gen long __tseq = _n
    quietly tsset __tseq
    newey `varlist', lag(`lag')
    restore
end

*==============================================================================*
* 1. CARGA + ASSERTS DE INTEGRIDAD (la corrida ABORTA si la base no es la canonica)
*==============================================================================*
import delimited "data/processed/base_cda.csv", clear varnames(1) case(preserve)

assert _N == 84
quietly count if missing(i_Gs)
assert r(N) == 0                              // primaria completa
quietly count if missing(i_USD)
assert r(N) == 0
quietly count if missing(i_Gs_pond)
assert r(N) == 1                              // ponderada: 1 faltante (2006Q4)
quietly count if missing(TC)
assert r(N) == 0
quietly count if missing(dolariz)
assert r(N) == 0
quietly count if missing(IPC_PY)
assert r(N) == 0
quietly count if missing(IPC_US)
assert r(N) == 0
quietly count if missing(tpm)
assert r(N) == 29                             // TPM nace 2011-05 (55/84)
di as result "[OK] Integridad de la base canonica verificada (84 filas; cobertura esperada por serie)"

gen t = yq(anio, trimestre)
format t %tq
tsset t

*==============================================================================*
* 2. VARIABLES (definiciones unicas del documento)
*==============================================================================*
gen difn = i_Gs - i_USD                       // diferencial nominal (pp anuales)
gen difq = difn/4                             // trimestralizado (lineal)
gen difq_comp = ((1+i_Gs/100)^0.25 - (1+i_USD/100)^0.25)*100   // conversion compuesta (M7)
gen dep     = (TC - L.TC)/L.TC * 100
gen depLead = (F.TC - TC)/TC * 100
gen DR     = difq - depLead
gen DRcons = difq - max(depLead, 0) if !missing(depLead)   // max(.,0)=0 en Stata: el if evita un valor espurio
gen infl_py  = (IPC_PY - L4.IPC_PY)/L4.IPC_PY * 100
gen infl_us  = (IPC_US - L4.IPC_US)/L4.IPC_US * 100
gen InflDiff = infl_py - infl_us
gen rp = difn - InflDiff                      // DIFERENCIAL REAL RESIDUAL (proxy de prima de riesgo)

* --- regimen de Metas de Inflacion: PRINCIPAL desde 2011T2 (18-may-2011) ---
gen D2011   = (anio > 2011) | (anio == 2011 & trimestre >= 2)
gen D2011t1 = (anio >= 2011)                  // sensibilidad (bloque 8.8)

* --- volatilidad cambiaria: sd movil 8 trimestres de dep (ddof=1, n>=2) ---
capture which rangestat
if _rc == 0 {
    rangestat (sd) VolTC = dep, interval(t -7 0)
}
else {
    gen VolTC = .
    forvalues i = 1/`=_N' {
        quietly summarize dep in `=max(1,`i'-7)'/`i'
        if r(N) >= 2 quietly replace VolTC = r(sd) in `i'
    }
}
gen Iliq = dolariz                            // dolarizacion/segmentacion (proxy; NO "iliquidez" a secas)

* --- centrado EN LA MUESTRA DEL MODELO (n=80: 2005Q1-2024Q4) ---
quietly summarize VolTC if !missing(rp, VolTC, Iliq)
gen VolTC_c = VolTC - r(mean)
quietly summarize Iliq if !missing(rp, VolTC, Iliq)
gen Iliq_c = Iliq - r(mean)
gen DxV = D2011*VolTC_c
gen DxI = D2011*Iliq_c

* --- series de sensibilidad ---
gen difn_fdt  = i_Gs_fdt - i_USD_fdt
gen difq_fdt  = difn_fdt/4
gen difn_pond = i_Gs_pond - i_USD_pond
gen rp_pond   = difn_pond - InflDiff

label var rp       "Diferencial real residual (proxy prima de riesgo, pp)"
label var DR       "Dif. de retorno trimestral (estandar)"
label var DRcons   "Dif. de retorno trimestral (conservador)"
label var InflDiff "Dif. de inflacion interanual (pp)"
label var VolTC    "Volatilidad cambiaria (sd movil 8T)"
label var Iliq     "Dolarizacion de depositos / segmentacion (%)"
label var D2011    "Regimen MdI (=1 desde 2011T2)"

*==============================================================================*
* BLOQUE 8.1 - Descriptivos (Tabla 8.1)
*==============================================================================*
di _n(2) "{hline 78}"
di " 8.1 Descriptivos de rp, DR y DRcons (corte 2011T2)"
di "{hline 78}"
tabstat rp DR DRcons, by(D2011) stat(n mean p50 sd min max skewness kurtosis) col(stat)
tabstat rp DR DRcons, stat(n mean p50 sd min max skewness kurtosis) col(stat)

tsline rp, lcolor(navy) lwidth(medthick) yline(0, lpattern(dash) lcolor(gs8))   ///
    xline(`=yq(2011,2)', lcolor(red) lwidth(medthick))                          ///
    title("Diferencial real residual rp") subtitle("marca roja = 2011T2 (MdI)") ///
    ytitle("pp anuales") xtitle("") name(g_rp, replace)
graph export "output/figuras/f1_rp_serie.png", replace width(1600)

tsline DR DRcons, lcolor(navy orange) lpattern(solid dash)                      ///
    yline(0, lpattern(dot) lcolor(gs8)) xline(`=yq(2011,2)', lcolor(red))       ///
    legend(order(1 "DR (estandar)" 2 "DRcons (conservador)"))                   ///
    title("DR y DRcons") name(g_dr, replace)
graph export "output/figuras/f2_DR_DRcons.png", replace width(1600)

*==============================================================================*
* BLOQUE 8.2 - Regresion de Fama / UIP (Tabla 8.2)
*   depLead = a + b*difq ; UIP <=> b=1. difq esta en la MISMA frecuencia
*   trimestral que depLead (la version anualizada implicaria b_UIP=0,25).
*   RESULTADO con la serie primaria: b=1,16 (total) y b=1 NO se rechaza en
*   ninguna muestra (R2 ~ 0): el diferencial NO predice la depreciacion.
*==============================================================================*
di _n(2) "{hline 78}"
di " 8.2 Fama/UIP (regresor trimestralizado difq; sensib.: compuesta y fdt)"
di "{hline 78}"
foreach rgr in difq difq_comp difq_fdt {
    di _n "--- regresor: `rgr' ---"
    newey_nogap depLead `rgr', lag($HAC)
    test `rgr' = 1
    newey_nogap depLead `rgr' if D2011==0, lag($HAC)
    test `rgr' = 1
    newey_nogap depLead `rgr' if D2011==1, lag($HAC)
    test `rgr' = 1
}

di _n "-- DRcons: H0 mu<=0 vs HA mu>0 (t HAC unilateral) --"
newey_nogap DRcons, lag($HAC)
scalar tt = _b[_cons]/_se[_cons]
di "   TOTAL:  media=" %7.4f _b[_cons] "  t=" %6.3f tt "  p_1cola(mu>0)=" %6.4f ttail(e(df_r), tt)
newey_nogap DRcons if D2011==0, lag($HAC)
scalar tt = _b[_cons]/_se[_cons]
di "   PRE:    media=" %7.4f _b[_cons] "  t=" %6.3f tt "  p_1cola(mu>0)=" %6.4f ttail(e(df_r), tt)
newey_nogap DRcons if D2011==1, lag($HAC)
scalar tt = _b[_cons]/_se[_cons]
di "   POST:   media=" %7.4f _b[_cons] "  t=" %6.3f tt "  p_1cola(mu>0)=" %6.4f ttail(e(df_r), tt)

*==============================================================================*
* BLOQUE 8.3 - Integracion y quiebre estructural (Tabla 8.3)
*   Diagnostico DECLARADO MIXTO: rp aparece ~I(1) en la ventana completa, pero
*   dentro de cada regimen no hay evidencia contra la estacionariedad (KPSS
*   p>0,10 en ambos; ADF post p<0,001): el aparente I(1) es el QUIEBRE DE NIVEL
*   de 2011 (Perron, 1989). No se afirma cointegracion formal (ver 8.7).
*==============================================================================*
di _n(2) "{hline 78}"
di " 8.3 Integracion (ADF + KPSS) y quiebre endogeno"
di "{hline 78}"
preserve
quietly keep if !missing(rp)
quietly gen long __tseq = _n
quietly tsset __tseq
foreach v in rp VolTC Iliq InflDiff DR {
    di _n "--- `v' (ventana completa) ---"
    capture noisily dfuller `v'
    capture noisily kpss `v', maxlag(4)
}
di _n "--- rp por SUBMUESTRA (argumento Perron 1989) ---"
capture noisily dfuller rp if D2011==0
capture noisily kpss rp if D2011==0, maxlag(4)
capture noisily dfuller rp if D2011==1
capture noisily kpss rp if D2011==1, maxlag(4)

di _n "-- Quiebre endogeno de rp: Bai-Perron (xtbreak) --"
capture noisily xtbreak test rp, breaks(3) hypothesis(3)
capture noisily xtbreak estimate rp, breaks(1)
di as text "   (resultado canonico: 1 quiebre en 2011Q2 - indice colapsado; mapear con list __tseq fecha)"
capture noisily xtbreak estimate rp, breaks(3)
restore

di _n "-- Diferencia de medias pre/post 2011T2 (HAC) --"
newey_nogap rp D2011, lag($HAC)
newey_nogap DR D2011, lag($HAC)

*==============================================================================*
* BLOQUE 8.4 - Descomposicion de Fisher (Tabla 8.4)
*   difn = InflDiff + rp. En NIVELES (pp anuales); NO en % porque rp<0 pre-2011.
*==============================================================================*
di _n(2) "{hline 78}"
di " 8.4 Fisher: difn = InflDiff + rp (niveles pre/post 2011T2)"
di "{hline 78}"
tabstat difn InflDiff rp if !missing(rp), by(D2011) stat(n mean sd) col(stat) longstub
di _n "-- t HAC de la diferencia pre/post --"
newey_nogap rp D2011, lag($HAC)
newey_nogap InflDiff D2011, lag($HAC)
newey_nogap difn D2011, lag($HAC)

*==============================================================================*
* BLOQUE 8.5 - MODELO CENTRAL (Tabla 8.5): rp centrado, D2011 desde T2
*   rp = b0 + b2 VolTC_c + b3 Iliq_c + b4 D2011 + b5 DxV + b6 DxI  (HAC lag 4)
*   VEREDICTO UNICO (los comentarios NO contradicen el Wald):
*     - NIVEL: b4 = +3,89 (p<0,001)  -> cambio de nivel fuerte.
*     - PENDIENTES: Wald b5=b6=0 -> F=2,00, p=0,142: NO se rechaza estabilidad.
*       El Chow completo es significativo POR EL INTERCEPTO, no por pendientes;
*       NO constituye evidencia de recomposicion interna.
*==============================================================================*
di _n(2) "{hline 78}"
di " 8.5 MODELO CENTRAL: rp ~ VolTC_c + Iliq_c + D2011(T2) + DxV + DxI"
di "{hline 78}"
newey_nogap rp VolTC_c Iliq_c D2011 DxV DxI, lag($HAC)
estimates store m_central
test DxV DxI                       // Wald de pendientes (contraste central)
test D2011 DxV DxI                 // Chow completo (nivel + pendientes)
scalar t_b5 = _b[DxV]/_se[DxV]
scalar t_b6 = _b[DxI]/_se[DxI]
di "   b5 (DxV) p_1cola(HA<0)=" %6.4f (1 - ttail(e(df_r), t_b5)) "   [conjetura: <0]"
di "   b6 (DxI) p_1cola(HA>0)=" %6.4f ttail(e(df_r), t_b6)       "   [conjetura: >0; el dato va en contra]"
capture noisily esttab m_central using "output/tablas/tabla_8_5_modelo_central.csv", ///
    replace se p nostar wide title("Tabla 8.5 - Modelo central (HAC lag 4)")

*==============================================================================*
* BLOQUE 8.6 - Dominancia estandarizada y composicion interna (Tabla 8.6)
*   Analisis DESCRIPTIVO (no contraste formal). Advertencia: el R2 del bloque
*   post-2011 es ~0,03: hay poco que repartir entre VolTC e Iliq.
*==============================================================================*
di _n(2) "{hline 78}"
di " 8.6 Dominancia estandarizada + composicion Shapley (descriptivo)"
di "{hline 78}"
egen zV = std(VolTC) if !missing(rp, VolTC, Iliq)
egen zI = std(Iliq)  if !missing(rp, VolTC, Iliq)
gen DzV = D2011*zV
gen DzI = D2011*zI
newey_nogap rp zV zI D2011 DzV DzI, lag($HAC)
lincom zV + DzV
lincom zI + DzI
lincom (zI + DzI) - (zV + DzV)

di _n "-- Shapley/LMG pre/post (2 regresores: formula cerrada) --"
foreach d in 0 1 {
    quietly regress rp VolTC if D2011==`d'
    scalar r2v = e(r2)
    quietly regress rp Iliq if D2011==`d'
    scalar r2i = e(r2)
    quietly regress rp VolTC Iliq if D2011==`d'
    scalar r2b = e(r2)
    scalar cv = 0.5*(r2v + r2b - r2i)
    scalar ci = 0.5*(r2i + r2b - r2v)
    di "   D2011=`d': n=" e(N) "  R2=" %5.3f r2b "  %VolTC=" %5.1f 100*cv/(cv+ci) "  %Iliq=" %5.1f 100*ci/(cv+ci)
}

* ventana movil de 20 trimestres (Figura f4)
preserve
quietly keep if !missing(rp, VolTC, Iliq)
gen pctVol = .
gen pctIliq = .
local W = 20
forvalues i = `W'/`=_N' {
    local j = `i' - `W' + 1
    quietly regress rp VolTC in `j'/`i'
    scalar r2v = e(r2)
    quietly regress rp Iliq in `j'/`i'
    scalar r2i = e(r2)
    quietly regress rp VolTC Iliq in `j'/`i'
    scalar r2b = e(r2)
    scalar cv = 0.5*(r2v + r2b - r2i)
    scalar ci = 0.5*(r2i + r2b - r2v)
    quietly replace pctVol  = 100*cv/(cv+ci) in `i'
    quietly replace pctIliq = 100*ci/(cv+ci) in `i'
}
gen t2 = yq(anio, trimestre)
format t2 %tq
quietly tsset t2
tsline pctVol pctIliq, lcolor(cranberry navy) yline(50, lpattern(dot))          ///
    xline(`=yq(2011,2)', lcolor(black) lpattern(dash))                          ///
    legend(order(1 "Riesgo cambiario (VolTC)" 2 "Dolarizacion/segmentacion"))   ///
    title("Composicion del diferencial real residual (ventana movil 20T)")     ///
    name(g_mov, replace)
graph export "output/figuras/f4_composicion_movil.png", replace width(1600)
restore

*==============================================================================*
* BLOQUE 8.7 - Diagnosticos (Tabla 8.7)
*   BG lags 3-4 y White significativos -> por eso TODA la inferencia es HAC.
*   Integracion mixta declarada (8.3); Engle-Granger canonico SIN dummies con
*   CV de MacKinnon NO rechaza no-cointegracion (t=-2,91 > cv5% -3,85): NO se
*   afirma cointegracion formal. La lectura del modelo descansa en (i) rp
*   estacionaria DENTRO de cada regimen y (ii) residuos del modelo con quiebre
*   estacionarios (ADF -5,0; referencia informal: los CV con quiebre son mas
*   exigentes). Limitacion declarada en el documento.
*==============================================================================*
di _n(2) "{hline 78}"
di " 8.7 Diagnosticos del OLS paralelo"
di "{hline 78}"
preserve
quietly keep if !missing(rp, VolTC_c, Iliq_c)
quietly gen long __tseq = _n
quietly tsset __tseq
regress rp VolTC_c Iliq_c D2011 DxV DxI
estat bgodfrey, lags(1 2 3 4)
estat imtest, white
estat ovtest
vif
predict resid_c, resid
quietly summarize resid_c, detail
scalar JB = r(N)/6*(r(skewness)^2 + (r(kurtosis)-3)^2/4)
di "   Jarque-Bera residuos = " %6.2f JB "  p = " %5.3f chi2tail(2, JB)
dfuller resid_c, lags(3)
di as text "   (ADF de residuos del modelo CON quiebre: referencia informal, ver nota del bloque)"
restore

*==============================================================================*
* BLOQUE 8.8 - Robustez y sensibilidad (Tabla 8.8)
*==============================================================================*
di _n(2) "{hline 78}"
di " 8.8 Robustez: T1, DRcons, L.Iliq, gradual, crisis, HAC lags, post, encaje, ponderada"
di "{hline 78}"

di _n "--- (a) D2011 desde T1 (sensibilidad de fecha; hallazgo M3) ---"
gen DxV_t1 = D2011t1*VolTC_c
gen DxI_t1 = D2011t1*Iliq_c
newey_nogap rp VolTC_c Iliq_c D2011t1 DxV_t1 DxI_t1, lag($HAC)
test DxV_t1 DxI_t1

di _n "--- (b) DRcons como dependiente ---"
newey_nogap DRcons VolTC_c Iliq_c D2011 DxV DxI, lag($HAC)
test DxV DxI

di _n "--- (c) dolarizacion rezagada (simultaneidad) ---"
gen L_Iliq = L.Iliq
quietly summarize L_Iliq if !missing(rp, VolTC, L_Iliq)
gen LIliq_c = L_Iliq - r(mean)
gen DxLI = D2011*LIliq_c
newey_nogap rp VolTC_c LIliq_c D2011 DxV DxLI, lag($HAC)
test DxV DxLI

di _n "--- (d) regimen gradual (rampa 2011T2-2016; 1 desde 2017) ---"
gen Reg = 0
replace Reg = ((anio-2011)*4 + trimestre - 2 + 1)/23 if (anio>2011 | (anio==2011 & trimestre>=2)) & anio<=2016
replace Reg = 1 if anio>=2017
replace Reg = min(Reg, 1)
gen RxV = Reg*VolTC_c
gen RxI = Reg*Iliq_c
newey_nogap rp VolTC_c Iliq_c Reg RxV RxI, lag($HAC)
test RxV RxI

di _n "--- (e) control por la crisis 2008Q2-2009Q3 ---"
gen Crisis0809 = (anio==2008 & trimestre>=2) | (anio==2009 & trimestre<=3)
newey_nogap rp VolTC_c Iliq_c D2011 DxV DxI Crisis0809, lag($HAC)
test DxV DxI

di _n "--- (f) sensibilidad al rezago HAC ---"
foreach L in 3 5 6 {
    di "   -- lag `L' --"
    newey_nogap rp VolTC_c Iliq_c D2011 DxV DxI, lag(`L')
    test DxV DxI
}

di _n "--- (g) submuestra post-2011 (sin dummies) ---"
newey_nogap rp VolTC Iliq if D2011==1, lag($HAC)

di _n "--- (h) encaje post-2011: brecha EFECTIVA ME-MN como proxy alternativa ---"
* La brecha de encaje es un fenomeno con cobertura verificada desde 2012Q4
* (parcial, declarada en la base). Se usa SOLO en la submuestra post-2011;
* fuente: docs/11_encaje_fuentes.md. UN SOLO BLOQUE, sin contradicciones.
capture confirm variable encaje_me_efec
if !_rc {
    gen brecha_encaje = encaje_me_efec - encaje_mn_efec
    quietly summarize brecha_encaje if D2011==1
    if r(N) >= 10 & r(sd) > 0 {
        newey_nogap rp VolTC brecha_encaje if D2011==1, lag($HAC)
    }
    else di as text "   [aviso] brecha sin variacion suficiente post-2011."
}
else di as text "   [aviso] columnas de encaje ausentes."

di _n "--- (i) SERIE PONDERADA (sensibilidad de definicion; hallazgo C1/C2) ---"
* Con la ponderada las pendientes SI cambian (Wald p<0,001) y el nivel se
* mantiene (+1,76; p<0,001). Se reporta como SENSIBILIDAD: tiene empalme
* 2010/2011 (~0,7-0,8pp definicional) y sesgo de composicion por plazos.
* El colapso de calendario de newey_nogap afecta a ESTA serie (hueco 2006Q4).
newey_nogap rp_pond VolTC_c Iliq_c D2011 DxV DxI, lag($HAC)
test DxV DxI
newey_nogap rp_pond D2011, lag($HAC)

*==============================================================================*
* BLOQUE 8.9 - TPM (descriptivo post-2011; Tabla 8.9)
*   La TPM nace con el regimen (may-2011): no entra como regresor de la ventana
*   completa (seria colineal con D2011). Fuente: BCP, TPM.xlsx oficial
*   (data/raw/bcp/TPM.xlsx; serie mensual -> promedio trimestral).
*==============================================================================*
di _n(2) "{hline 78}"
di " 8.9 TPM y tasa CDA en guaranies (descriptivo, post-2011)"
di "{hline 78}"
quietly corr i_Gs tpm if !missing(tpm)
di "   corr(i_Gs, TPM) = " %6.3f r(rho) "   (n=" r(N) ")"
gen spread_GsTPM = i_Gs - tpm
summarize spread_GsTPM if !missing(tpm)
tsline i_Gs tpm if !missing(tpm), lcolor(navy cranberry) lpattern(solid dash)   ///
    legend(order(1 "CDA <=365d en Gs" 2 "TPM")) title("Tasa CDA en Gs y TPM")   ///
    name(g_tpm, replace)
graph export "output/figuras/f5_tpm_iGs.png", replace width(1600)

*==============================================================================*
* CIERRE
*==============================================================================*
di _n(2) "{hline 78}"
di " FIN DE LA CORRIDA CANONICA."
di " Log:     output/log_tesis_cda.log"
di " Tablas:  output/tablas/   |  Figuras: output/figuras/"
di " Validacion cruzada independiente (Python): python/corrida_canonica.py"
di "   -> output/canonica/resultados.json  (misma especificacion, mismos numeros)"
di "{hline 78}"
log close
