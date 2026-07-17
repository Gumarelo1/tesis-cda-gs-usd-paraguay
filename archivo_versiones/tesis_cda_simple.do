*==============================================================================*
*  TESIS: Diferencial de Retornos CDA Guaranies vs Dolares (Paraguay 2004-2024)
*  Scavone & Fernandez - VERSION SIMPLE PARA LA DEFENSA
*------------------------------------------------------------------------------*
*  Este do-file produce los resultados centrales de la tesis con codigo
*  directo, pensado para poder explicarse linea por linea.
*  La version completa (tesis_cda.do) agrega robustez extra: series de fin de
*  trimestre, descomposicion de Shapley y evolucion por ventanas moviles.
*  En todo lo que comparten, AMBOS do-files dan exactamente los mismos numeros.
*
*  COMO CORRER: 1) ajustar la ruta de abajo si hace falta; 2) do este archivo.
*==============================================================================*

clear all
set more off
version 17

cd "C:/Users/Laviero Scavone/Desktop/tes"          // <- unica edicion manual
capture mkdir "output"
capture mkdir "output/figuras"

* Paquetes (solo la primera vez; requieren internet)
capture which rangestat
if _rc capture ssc install rangestat, replace
capture which xtbreak
if _rc capture ssc install xtbreak, replace
capture which kpss
if _rc capture ssc install kpss, replace

*==============================================================================*
* 1. CARGAR LA BASE Y CONSTRUIR LAS VARIABLES
*==============================================================================*
import delimited "data/processed/base_cda.csv", clear varnames(1) case(preserve)

gen t = yq(anio, trimestre)        // fecha trimestral
format t %tq
tsset t

* Diferencial de tasas (puntos porcentuales)
gen difn = i_Gs - i_USD            // anual
gen difq = difn/4                  // trimestralizado (misma frecuencia que la depreciacion)

* Depreciacion del guarani (%): de este trimestre y del siguiente
gen dep     = (TC - L.TC)/L.TC * 100
gen depLead = (F.TC - TC)/TC * 100

* Series de retorno de la tesis (ecuaciones 5 y 6)
gen DR     = difq - depLead                                  // estandar
gen DRcons = difq - max(depLead, 0) if !missing(depLead)     // conservadora

* Regimen de Metas de Inflacion (desde 2011T1)
gen D2011 = (anio >= 2011)

* Diferencial de inflacion interanual (Fisher) y prima de riesgo (ecuacion 8)
gen infl_py  = (IPC_PY - L4.IPC_PY)/L4.IPC_PY * 100
gen infl_us  = (IPC_US - L4.IPC_US)/L4.IPC_US * 100
gen InflDiff = infl_py - infl_us
gen rp = difn - InflDiff           // prima de riesgo = diferencial - inflacion

* Volatilidad cambiaria: desvio movil de 8 trimestres de la depreciacion
rangestat (sd) VolTC = dep, interval(t -7 0)

* Iliquidez (proxy: dolarizacion de depositos) e interacciones con el regimen
gen Iliq    = dolariz
gen D_VolTC = D2011 * VolTC
gen D_Iliq  = D2011 * Iliq

* Rezago de la dolarizacion construido AQUI, sobre el calendario completo (antes
* del drop), para que no cruce los huecos de CDA. Se usa en la robustez de
* simultaneidad (Etapa 5-TER (5)). La dolarizacion no tiene huecos, asi que
* L.Iliq solo falta en 2004T1.
gen L_Iliq = L.Iliq

*------------------------------------------------------------------------------*
* MANEJO DE DATOS FALTANTES (importante para la defensa):
* Con la serie primaria del tramo CDA <=365 dias (fuente unica, sin empalme)
* NO quedan huecos trimestrales: el 'drop' de abajo no elimina filas y se
* conserva solo como salvaguarda. Los faltantes de borde (inicio 2004 por
* rezagos, fin 2024 por el adelanto) no molestan a 'newey'. La serie de
* sensibilidad i_*_pond (promedio ponderado) si tiene un hueco (2006T4); su
* bloque al final lo maneja con preserve/drop/renumeracion.
*------------------------------------------------------------------------------*
drop if missing(i_Gs)              // 0 filas con la serie <=365 (sin huecos); salvaguarda
gen tt = _n
tsset tt

*==============================================================================*
* ETAPA 1 - Estadistica descriptiva de las series de retorno
*   Que muestra: DR promedio cae de 1.70 (pre-2011) a -0.32 (post-2011);
*   la serie conservadora es negativa en todos los cortes.
*==============================================================================*
di _n "======= ETAPA 1: descriptivas ======="
tabstat DR DRcons, by(D2011) stat(n mean p50 sd min max) col(stat)

tsline DR, yline(0, lpattern(dash)) ///
    title("Diferencial de retorno DR (Gs vs USD)") ytitle("pp trimestral") ///
    name(gDR, replace)
graph export "output/figuras/simple_DR.png", replace width(1200)

*==============================================================================*
* ETAPA 2 - Regresion de Fama + prueba de la serie conservadora
*   Que muestra: beta positivo pero impreciso; beta=1 (UIP) NO se rechaza en
*   ninguna muestra (total p=0.90; R2~0.01 -> contraste poco informativo).
*   La media de DRcons no es positiva.
*==============================================================================*
di _n "======= ETAPA 2: Fama (UIP <=> beta=1) ======="
newey depLead difq, lag(4)
test difq = 1                       // total: p=0.898 -> beta=1 NO se rechaza (baja potencia)
newey depLead difq if D2011==0, lag(4)
test difq = 1
newey depLead difq if D2011==1, lag(4)
test difq = 1

di _n "--- Media de DRcons (H0: <=0). El coeficiente _cons ES la media ---"
newey DRcons, lag(4)                // total:    -1.00 (no positiva)
di "  p 1-cola (mu>0) = " %6.4f ttail(e(df_r), _b[_cons]/_se[_cons])
newey DRcons if D2011==0, lag(4)    // pre-2011: -0.61
di "  p 1-cola (mu>0) = " %6.4f ttail(e(df_r), _b[_cons]/_se[_cons])
newey DRcons if D2011==1, lag(4)    // post:     -1.20 (significativamente < 0)
di "  p 1-cola (mu>0) = " %6.4f ttail(e(df_r), _b[_cons]/_se[_cons])
* porcentaje de trimestres con DRcons > 0 (lo reporta el documento)
count if DRcons>0 & !missing(DRcons)
count if !missing(DRcons)

*==============================================================================*
* ETAPA 3 - Estacionariedad y quiebre estructural
*   Que muestra: DR es estacionaria (sin riesgo de regresion espuria) y el
*   quiebre endogeno cae en 2008T2 (crisis global), no en 2011.
*==============================================================================*
di _n "======= ETAPA 3: raiz unitaria y quiebre ======="
quietly varsoc DR, maxlag(8)        // seleccion de rezagos por BIC (elige 0)
dfuller DR, lags(0)                 // ADF (nula=raiz unitaria): -7.42 -> estacionaria
capture noisily kpss DR, maxlag(4)  // KPSS (nula=estacionaria, invierte la nula): confirma
pperron DR                          // PP como referencia adicional: -7.62

newey DR D2011, lag(4)              // dif. de medias pre/post: -2.02 (p=0.058)

quietly regress DR L.DR
capture noisily estat sbsingle      // quiebre unico endogeno -> 2008T2
capture noisily estat sbcusum
capture noisily xtbreak test DR, breaks(3) hypothesis(3)   // UDmax/WDmax
capture noisily xtbreak estimate DR, breaks(3)              // Bai-Perron -> 2008T2 (+2009T3, 2011T3)

*==============================================================================*
* ETAPA 4 - Descomposicion de Fisher (CAPA 1 de la hipotesis: CONFIRMADA)
*   Que muestra: pre-2011 la prima de riesgo era NEGATIVA (-2.18 pp: el
*   diferencial ni compensaba la inflacion); tras 2011 emerge una prima
*   positiva (+3.70, p<0.001), la prima de inflacion cae (-2.83, p<0.001)
*   y el diferencial nominal sube (+0.87, p=0.018).
*==============================================================================*
di _n "======= ETAPA 4: Fisher (inflacion vs riesgo) ======="
tabstat difn InflDiff rp, by(D2011) stat(n mean sd) col(stat)

* niveles (pp anuales) dentro del diferencial. NO se reporta en % porque la
* prima de riesgo es NEGATIVA pre-2011 (-2.18): una descomposicion porcentual
* daria valores sin sentido (p.ej. 190% inflacion / -90% riesgo).
foreach d in 0 1 {
    quietly summarize difn if D2011==`d' & !missing(rp)
    local m_difn = r(mean)
    quietly summarize InflDiff if D2011==`d' & !missing(rp)
    local m_infl = r(mean)
    quietly summarize rp if D2011==`d' & !missing(rp)
    local m_rp = r(mean)
    di "  D2011=`d':  difn=" %5.2f `m_difn' "  prima inflacion=" %5.2f `m_infl' ///
       "  prima riesgo=" %5.2f `m_rp' "  (pp anuales)"
}

* significancia del cambio (el coef. de D2011 es la diferencia post-pre)
newey InflDiff D2011, lag(4)        // -2.83 (p=0.0005): la prima de inflacion CAE
newey rp       D2011, lag(4)        // +3.70 (p<0.001): la prima SUBE (de negativa a positiva)

*==============================================================================*
* ETAPA 5 - Modelo con interacciones sobre la prima de riesgo rp
*   NOTA METODOLOGICA (circularidad de Fisher): por construccion el diferencial
*   contiene al diferencial de inflacion (difn = InflDiff + rp). Por eso NO se
*   regresa DR (ni difn) usando InflDiff como regresor: eso dejaria la inflacion
*   en AMBOS lados de la ecuacion. El contraste de determinantes se hace SOBRE
*   rp, donde la inflacion ya fue removida del lado derecho. (Ajuste central
*   pedido por el tutor.)
*   Que muestra: los coeficientes direccionales b5 y b6 no son significativos
*   -> la conjetura especifica (migracion hacia la iliquidez) NO se sostiene.
*   La colinealidad D2011~D_Iliq (corr 0.99) limita los contrastes individuales;
*   el contraste conjunto (general) esta en la Etapa 5-BIS.
*==============================================================================*
di _n "======= ETAPA 5: interacciones sobre la prima de riesgo rp ======="
newey rp VolTC Iliq D2011 D_VolTC D_Iliq, lag(4)
test D_VolTC                        // b5 (esperado <0): -0.19, p=0.48
di "  p 1-cola (b5<0) = " %6.4f 1-ttail(e(df_r), _b[D_VolTC]/_se[D_VolTC])
test D_Iliq                         // b6 (esperado >0): -0.13, p=0.43 (n.s.)
di "  p 1-cola (b6>0) = " %6.4f ttail(e(df_r), _b[D_Iliq]/_se[D_Iliq])

newey DRcons VolTC Iliq D2011 D_VolTC D_Iliq, lag(4)   // robustez (sin InflDiff)

*==============================================================================*
* ETAPA 5-BIS - CONTRASTE CENTRAL DE LA TESIS (forma general, CONFIRMADA)
*   Hipotesis: "tras 2011 la composicion de la prima de riesgo CAMBIO"
*   (sin imponer direccion). Se contrasta sobre rp con una prueba de Wald
*   CONJUNTA de las interacciones, robusta a la colinealidad que anula los
*   contrastes individuales.
*   Que muestra: Wald F=0.36 (p=0.698) -> las PENDIENTES no cambian; Chow
*   F=16.81 (p<0.001) -> el cambio es de NIVEL (salto de la prima). La
*   recomposicion ocurre entre primas (Fisher), no dentro de la prima.
*==============================================================================*
di _n "======= ETAPA 5-BIS: cambio de composicion de la prima de riesgo ======="
newey rp VolTC Iliq D2011 D_VolTC D_Iliq, lag(4)
test D_VolTC D_Iliq                 // pendientes: F=0.36, p=0.698 (NO cambian)
test D2011 D_VolTC D_Iliq           // Chow completo: F=16.81, p<0.001 (nivel)

di _n "--- La composicion antes y despues (coeficientes por subperiodo) ---"
newey rp VolTC Iliq if D2011==0, lag(4)   // pre : VolTC=+0.27 (p=0.017); Iliq n.s.
newey rp VolTC Iliq if D2011==1, lag(4)   // post: ambos n.s. (R2~0)

*==============================================================================*
* ETAPA 5-TER - AJUSTES DE LA REVISION METODOLOGICA (devolucion del tutor)
*   (1) Modelo principal: dependiente = prima de riesgo rp (sin InflDiff en el
*       lado derecho: elimina la circularidad de Fisher) y variables CENTRADAS
*       en su media (asi D2011 se interpreta en valores promedio).
*   (2) Prueba de DOMINANCIA post-2011 en unidades comparables (estandarizada):
*       el contraste directo de la "migracion" que promete el titulo.
*   (3) Orden de integracion de TODAS las series (ADF + KPSS).
*   (4) Normalidad de residuos (Jarque-Bera).
*   (5) Robustez: dolarizacion rezagada (mitiga la simultaneidad declarada).
*==============================================================================*
di _n "======= ETAPA 5-TER (1): modelo principal revisado (rp, centrado) ======="
capture which kpss
if _rc capture ssc install kpss, replace

quietly summarize VolTC if !missing(rp, VolTC, Iliq)
gen VolTC_c = VolTC - r(mean)
quietly summarize Iliq if !missing(rp, VolTC, Iliq)
gen Iliq_c = Iliq - r(mean)
gen DxV = D2011*VolTC_c
gen DxI = D2011*Iliq_c
newey rp VolTC_c Iliq_c D2011 DxV DxI, lag(4)   // D2011 = +3.71 (p<0.001): SALTO DE NIVEL
test DxV DxI          // Wald conjunto: F=0.36, p=0.698 (pendientes estables)
test D2011 DxV DxI    // Chow completo:                 F=16.81, p<0.001

di _n "======= ETAPA 5-TER (2): dominancia post-2011 (estandarizada) ======="
egen zV = std(VolTC) if !missing(rp, VolTC, Iliq)
egen zI = std(Iliq)  if !missing(rp, VolTC, Iliq)
gen DzV = D2011*zV
gen DzI = D2011*zI
newey rp zV zI D2011 DzV DzI, lag(4)
lincom zV + DzV                    // incidencia post riesgo cambiario: +0.20 sd
lincom zI + DzI                    // incidencia post iliquidez:        +0.08 sd
lincom (zI + DzI) - (zV + DzV)     // DOMINANCIA: t=-0.29, p=0.772 ->
*   incidencias indistinguibles post-2011 (sin migracion de dominancia)

di _n "======= ETAPA 5-TER (3): integracion de todas las series (ADF + KPSS) ======="
foreach v in rp VolTC Iliq InflDiff {
    di _n "--- `v' ---"
    dfuller `v' if !missing(rp)
    capture noisily kpss `v' if !missing(rp), maxlag(4)
}
*   Hallazgo: las explicativas se comportan como no estacionarias (la
*   dolarizacion cae de forma sostenida) -> limite declarado; ver (5).

di _n "======= ETAPA 5-TER (4): Jarque-Bera de los residuos ======="
quietly regress rp VolTC_c Iliq_c D2011 DxV DxI
predict resid_rp, resid
quietly summarize resid_rp, detail
scalar JB = r(N)/6*(r(skewness)^2 + (r(kurtosis)-3)^2/4)
di "  JB = " %6.2f JB "   p = " %5.3f chi2tail(2, JB) "   (menor a 5.99 => normalidad no rechazada)"
* cointegracion de Engle-Granger: ADF (3 rezagos) sobre los residuos = -4.77
* (p<0.001) -> residuos estacionarios -> COINTEGRAN -> la relacion NO es espuria
dfuller resid_rp, lags(3)

di _n "======= ETAPA 5-TER (5): robustez con dolarizacion rezagada ======="
* L_Iliq ya se construyo arriba sobre el calendario (no cruza los huecos)
quietly summarize L_Iliq if !missing(rp, VolTC, L_Iliq)
gen LIliq_c = L_Iliq - r(mean)
gen DxLI = D2011*LIliq_c
newey rp VolTC_c LIliq_c D2011 DxV DxLI, lag(4)
test DxV DxLI         // mismo diagnostico: F=0.40, p=0.673 (estables)

di _n "======= ETAPA 5-TER (6): robustez al supuesto de regimen 0/1 (gradual) ======="
* El profe observa que la credibilidad del regimen se consolido de a poco.
* Se re-testea con un indice de INTENSIDAD del regimen que sube gradualmente
* (0 pre-2011; rampa lineal 2011-2016; 1 desde 2017) en lugar de la dummy 0/1.
* Si la recomposicion se sostiene, no depende del supuesto de salto instantaneo.
gen Reg = 0
replace Reg = ((anio-2011)*4 + (trimestre-1) + 1)/24 if anio>=2011 & anio<=2016
replace Reg = 1 if anio>=2017
gen RxV = Reg*VolTC_c
gen RxI = Reg*Iliq_c
newey rp VolTC_c Iliq_c Reg RxV RxI, lag(4)
test RxV RxI          // gradual: F=2.26, p=0.111 (conjunta n.s.; RxI -0.49, p=0.038)

*==============================================================================*
* ETAPA 6 - Diagnosticos y robustez
*   Que muestra: hay heterocedasticidad y algo de autocorrelacion -> por eso
*   todo se estima con errores HAC (Newey-West). La forma funcional esta bien.
*   El VIF confirma la colinealidad D2011~D_Iliq (limitacion declarada).
*==============================================================================*
di _n "======= ETAPA 6: diagnosticos ======="
regress rp VolTC Iliq D2011 D_VolTC D_Iliq
estat bgodfrey, lags(1 2 3 4)       // autocorrelacion (lag 3-4 significativa)
estat imtest, white                 // heterocedasticidad (p<0.001)
estat ovtest                        // RESET: forma funcional OK (p=0.88)
vif                                 // D2011 y D_Iliq ~180-200: colinealidad
correlate D2011 D_Iliq D_VolTC Iliq

* Robustez: submuestra post-2011 (sin dummies)
newey rp     VolTC Iliq if D2011==1, lag(4)
newey DRcons VolTC Iliq if D2011==1, lag(4)

* Nota: la robustez con brecha de encaje NO es factible (la serie oficial
* homogenea 2004-2024 no existe en las publicaciones publicas del BCP).

*==============================================================================*
* ETAPA 7 (SENSIBILIDAD) - Serie alternativa: promedio ponderado de plazos
*   La serie primaria es el tramo CDA <=365 dias (fuente unica, sin empalme).
*   La ponderada (i_*_pond: IF 2004-2010 + Cuadro 31 2011-2024, con empalme y
*   sesgo de composicion por plazos) se reestima SOLO como sensibilidad:
*   el salto de NIVEL se mantiene (+1.80, p=0.021) pero las interacciones
*   aparecen significativas (Wald F=6.10, p=0.004). Conclusion declarada:
*   el cambio de nivel es robusto a la definicion de la serie; la evidencia
*   de pendientes es sensible a ella.
*==============================================================================*
di _n "======= ETAPA 7: sensibilidad a la definicion de la serie ======="
preserve
gen difn_pond = i_Gs_pond - i_USD_pond
gen rp_pond   = difn_pond - InflDiff
drop if missing(rp_pond, VolTC, Iliq)     // la ponderada tiene hueco en 2006T4
gen tt2 = _n
tsset tt2
newey rp_pond VolTC Iliq D2011 D_VolTC D_Iliq, lag(4)
test D_VolTC D_Iliq   // ponderada: F=6.10, p=0.004 (vs 0.36/0.698 con <=365)
newey rp_pond D2011, lag(4)               // salto de nivel ponderada: +1.80 (p=0.021)
restore
tsset tt

di _n "======= FIN. Version completa con robustez extra: tesis_cda.do ======="
