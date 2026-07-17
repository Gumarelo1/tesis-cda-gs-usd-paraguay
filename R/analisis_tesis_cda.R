# ============================================================================
#  TESIS CDA Gs/USD 2004-2024  —  ANALISIS ECONOMETRICO EN R
#  Reproduccion independiente del do-file de Stata (stata/tesis_cda.do).
#  Serie principal: tramo CDA <=365 dias. Corte institucional principal: 2011:T2.
#  Autores: Laviero Scavone y Matias Fernandez.
#
#  Objetivo: demostrar, con datos reales y sin ajustes, los resultados de la
#  tesis, y verificar que R reproduce las cifras de Stata (control automatico).
#
#  COMO CORRERLO (VS Code, extension R + radian):
#    1) Abrir esta carpeta del proyecto como workspace.
#    2) Ejecutar el script entero (o por bloques con Ctrl+Enter).
#    Lee: data/processed/base_cda.csv   Escribe: output_R/
#
#  CONVENCION HAC: se replica el 'newey' de Stata con
#    sandwich::NeweyWest(., lag=4, prewhite=FALSE, adjust=TRUE)  [Bartlett + n/(n-k)]
#  y p-valores con distribucion t(n-k) (coeftest con df=residual).
# ============================================================================

## ---- 0. Paquetes -----------------------------------------------------------
pkgs <- c("sandwich","lmtest","car","tseries","strucchange","zoo","ggplot2","writexl")
faltan <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(faltan)) install.packages(faltan, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, require, character.only = TRUE))

RUTA_BASE <- "data/processed/base_cda.csv"     # ajustar si se corre desde otra carpeta
if (!file.exists(RUTA_BASE)) stop("No encuentro ", RUTA_BASE, " . Abri la carpeta del proyecto como workspace.")
dir.create("output_R", showWarnings = FALSE)
dir.create("output_R/figuras", showWarnings = FALSE)
HAC <- 4L

## Helper: regresion con errores HAC a la Stata (Newey-West) --------------------
newey_stata <- function(formula, data, lag = HAC) {
  m  <- lm(formula, data = data)                       # na.omit descarta filas con NA (listwise)
  V  <- sandwich::NeweyWest(m, lag = lag, prewhite = FALSE, adjust = TRUE)
  ct <- lmtest::coeftest(m, vcov. = V, df = df.residual(m))
  list(m = m, V = V, ct = ct, n = nobs(m), dfr = df.residual(m))
}
waldF <- function(fit, terms) {
  h <- paste0(terms, " = 0")
  lh <- car::linearHypothesis(fit$m, h, vcov. = fit$V, test = "F")
  c(F = lh$F[2], p = lh$`Pr(>F)`[2], df1 = length(terms), df2 = fit$dfr)
}
p2  <- function(b, se, dfr) 2 * pt(abs(b/se), dfr, lower.tail = FALSE)  # p bilateral t

## ---- 1. Carga e INTEGRIDAD -------------------------------------------------
d <- read.csv(RUTA_BASE, check.names = FALSE, stringsAsFactors = FALSE)
d <- d[order(d$anio, d$trimestre), ]
stopifnot(nrow(d) == 84)                                            # 84 trimestres
stopifnot(all(diff(d$anio * 4 + d$trimestre) == 1))                # contiguos
stopifnot(sum(is.na(d$i_Gs)) == 0, sum(is.na(d$i_USD)) == 0)       # serie principal completa
cat("[OK] Integridad: 84 obs contiguas, serie principal sin faltantes.\n")

## ---- 2. Construccion de variables (identica al do-file, corte 2011:T2) ------
lagn  <- function(x, k = 1) c(rep(NA, k), head(x, -k))
leadn <- function(x, k = 1) c(tail(x, -k), rep(NA, k))
dplyr_lag <- function(x, k = 1) lagn(x, k)

d$dep     <- 100 * (d$TC / lagn(d$TC) - 1)
d$depLead <- 100 * (leadn(d$TC) / d$TC - 1)
d$difn    <- d$i_Gs - d$i_USD
d$difq    <- d$difn / 4
d$DR      <- d$difq - d$depLead
d$DRcons  <- ifelse(is.na(d$depLead), NA, d$difq - pmax(d$depLead, 0))
d$infl_py <- 100 * (d$IPC_PY / lagn(d$IPC_PY, 4) - 1)
d$infl_us <- 100 * (d$IPC_US / lagn(d$IPC_US, 4) - 1)
d$InflDiff<- d$infl_py - d$infl_us
d$rp      <- d$difn - d$InflDiff
d$D2011   <- as.integer(d$anio > 2011 | (d$anio == 2011 & d$trimestre >= 2))   # 2011:T2
d$D2011_T1<- as.integer(d$anio >= 2011)
d$crisis  <- as.integer((d$anio == 2008 & d$trimestre >= 2) | (d$anio == 2009 & d$trimestre <= 3))
# volatilidad cambiaria: sd movil de 8 trimestres, min 2, ddof=1 (== rangestat)
volsd <- function(x){ x <- x[!is.na(x)]; if (length(x) >= 2) sd(x) else NA_real_ }
d$VolTC <- zoo::rollapplyr(d$dep, width = 8, FUN = volsd, partial = TRUE, fill = NA)
d$dolar <- d$dolariz
# chequeo de identidad de Fisher y de formula de inflacion
stopifnot(all(abs(d$infl_py - d$infl_PY_yoy) < 1e-3, na.rm = TRUE))
stopifnot(all(abs(d$difn - (d$InflDiff + d$rp)) < 1e-6, na.rm = TRUE))
# centrado sobre la muestra del modelo
sm <- with(d, !is.na(rp) & !is.na(VolTC) & !is.na(dolar))
d$VolTC_c <- d$VolTC - mean(d$VolTC[sm]); d$dolar_c <- d$dolar - mean(d$dolar[sm])
d$DxVolTC <- d$D2011 * d$VolTC_c; d$DxDolar <- d$D2011 * d$dolar_c
# variantes fin de trimestre y ponderada
d$difn_fdt <- d$i_Gs_fdt - d$i_USD_fdt; d$rp_fdt <- d$difn_fdt - d$InflDiff
d$difn_pond<- d$i_Gs_pond - d$i_USD_pond; d$rp_pond <- d$difn_pond - d$InflDiff
cat("[OK] Variables construidas (Fisher identity y formula de inflacion verificadas).\n\n")

RES <- list()   # acumulador de resultados para exportar

## ---- 3. Descriptivos por subperiodo ---------------------------------------
descr <- do.call(rbind, lapply(c("difn","DR","DRcons","InflDiff","rp"), function(v){
  do.call(rbind, lapply(list(pre = d$D2011==0, post = d$D2011==1, total = rep(TRUE,nrow(d))), function(f){
    x <- d[[v]][f]; data.frame(var=v, n=sum(!is.na(x)), media=mean(x,na.rm=TRUE),
      sd=sd(x,na.rm=TRUE), min=min(x,na.rm=TRUE), p50=median(x,na.rm=TRUE), max=max(x,na.rm=TRUE))
  }))
}))
RES$descriptivos <- descr
cat("== Descriptivos (media por subperiodo) ==\n")
print(aggregate(cbind(difn,DR,DRcons,InflDiff,rp) ~ D2011, d, mean, na.rm=TRUE), digits=4)

## ---- 4. Regresion de Fama (UIP: beta=1) ------------------------------------
fama <- do.call(rbind, lapply(list(Total=rep(TRUE,nrow(d)), Pre=d$D2011==0, Post=d$D2011==1), function(f){
  f2 <- newey_stata(depLead ~ difq, data = d[f,])
  b <- coef(f2$m)["difq"]; se <- sqrt(f2$V["difq","difq"])
  lh <- car::linearHypothesis(f2$m, "difq = 1", vcov.=f2$V, test="F")
  data.frame(n=f2$n, beta=b, se=se, p_b0=p2(b,se,f2$dfr), p_b1=lh$`Pr(>F)`[2])
}))
rownames(fama) <- c("Total","Pre","Post"); RES$fama <- fama
cat("\n== Fama: beta y p(beta=1) ==\n"); print(fama, digits=4)

## ---- 5. Prueba conservadora DRcons -----------------------------------------
drc <- do.call(rbind, lapply(list(Total=rep(TRUE,nrow(d)), Post=d$D2011==1), function(f){
  f2 <- newey_stata(DRcons ~ 1, data=d[f,]); b <- coef(f2$m)[1]; se <- sqrt(f2$V[1,1])
  data.frame(n=f2$n, media=b, p_1cola_H1mayor0 = pt(b/se, f2$dfr, lower.tail=FALSE))
}))
rownames(drc) <- c("Total","Post"); RES$DRcons <- drc
cat("\n== DRcons (H1: media>0) ==\n"); print(drc, digits=4)

## ---- 6. Descomposicion de Fisher -------------------------------------------
fis <- do.call(rbind, lapply(c("rp","InflDiff","difn"), function(v){
  f2 <- newey_stata(as.formula(paste(v,"~ D2011")), data=d[!is.na(d$rp),])
  b <- coef(f2$m)["D2011"]; se <- sqrt(f2$V["D2011","D2011"])
  data.frame(componente=v, delta=b, p=p2(b,se,f2$dfr))
}))
RES$fisher <- fis
cat("\n== Fisher: cambio de medias pre/post 2011:T2 ==\n"); print(fis, digits=4)

## ---- 7. MODELO CENTRAL MCO-HAC (nucleo del contraste) ----------------------
m0 <- newey_stata(rp ~ VolTC_c + dolar_c + D2011 + DxVolTC + DxDolar, data = d)
cat("\n== MODELO CENTRAL (rp, corte 2011:T2) ==\n"); print(m0$ct)
w_slopes <- waldF(m0, c("DxVolTC","DxDolar"))
w_chow   <- waldF(m0, c("D2011","DxVolTC","DxDolar"))
cat(sprintf("Wald pendientes: F=%.3f  p=%.3f   |   Chow: F=%.2f  p=%.3f\n",
            w_slopes["F"], w_slopes["p"], w_chow["F"], w_chow["p"]))
RES$modelo_central <- as.data.frame(m0$ct[]); RES$modelo_central$termino <- rownames(m0$ct[])

## ---- 8. Diagnosticos (sobre el OLS clasico) --------------------------------
ols <- lm(rp ~ VolTC_c + dolar_c + D2011 + DxVolTC + DxDolar, data = d)
bg  <- sapply(1:4, function(L) lmtest::bgtest(ols, order=L)$p.value)
wh  <- lmtest::bptest(ols, ~ fitted(ols) + I(fitted(ols)^2))       # White (version fitted)
rst <- lmtest::resettest(ols, power=2:3, type="fitted")
jb  <- tseries::jarque.bera.test(residuals(ols))
vif <- car::vif(ols)
cat("\n== Diagnosticos ==\n")
cat(sprintf("Breusch-Godfrey p (lags 1-4): %s\n", paste(round(bg,4), collapse=", ")))
cat(sprintf("White p=%.4f | RESET p=%.4f | Jarque-Bera p=%.4f\n", wh$p.value, rst$p.value, jb$p.value))
cat("VIF:\n"); print(round(vif,2))

## ---- 9. Raices unitarias y quiebre estructural -----------------------------
adf_kpss <- function(v){
  x <- na.omit(d[[v]])
  adf <- tryCatch(suppressWarnings(tseries::adf.test(x, k=4)$p.value), error=function(e) NA)
  kps <- tryCatch(suppressWarnings(tseries::kpss.test(x)$p.value), error=function(e) NA)
  data.frame(var=v, ADF_p=adf, KPSS_p=kps)
}
RES$raiz_unitaria <- do.call(rbind, lapply(c("DR","rp","VolTC","dolar","InflDiff"), adf_kpss))
cat("\n== Raiz unitaria (ADF: H0 raiz unit.; KPSS: H0 estacionaria) ==\n"); print(RES$raiz_unitaria, digits=3)
# quiebre endogeno de la prima rp (sup-Wald tipo Quandt-Andrews)
rp_ts <- na.omit(d$rp); fs <- strucchange::Fstats(rp_ts ~ 1, from=0.15)
bp <- which.max(fs$Fstats); fecha_bp <- (which(!is.na(d$rp))[1] + bp - 1)  # posicion aprox
cat(sprintf("Quiebre endogeno de rp: sup-F en la observacion ~%d (esperado 2011:T2); p=%.4f\n",
            bp, strucchange::sctest(fs)$p.value))

## ---- 10. Robustez ----------------------------------------------------------
rob_row <- function(nombre, formula, data, terms, levelvar){
  f2 <- newey_stata(formula, data=data); b <- coef(f2$m)[levelvar]; se <- sqrt(f2$V[levelvar,levelvar])
  w  <- waldF(f2, terms)
  data.frame(especificacion=nombre, n=f2$n, nivel=b, p_nivel=p2(b,se,f2$dfr), Wald_F=w["F"], Wald_p=w["p"])
}
# T1
d$VolTC_c1 <- d$VolTC - mean(d$VolTC[sm]); d$dolar_c1 <- d$dolar - mean(d$dolar[sm])
d$DxV_T1 <- d$D2011_T1*d$VolTC_c1; d$DxD_T1 <- d$D2011_T1*d$dolar_c1
rob <- rbind(
  rob_row("Principal 2011:T2", rp ~ VolTC_c+dolar_c+D2011+DxVolTC+DxDolar, d, c("DxVolTC","DxDolar"), "D2011"),
  rob_row("Sensib. 2011:T1", rp ~ VolTC_c1+dolar_c1+D2011_T1+DxV_T1+DxD_T1, d, c("DxV_T1","DxD_T1"), "D2011_T1")
)
# fin de trimestre
smf <- with(d, !is.na(rp_fdt) & !is.na(VolTC) & !is.na(dolar))
d$VolTC_cf <- d$VolTC-mean(d$VolTC[smf]); d$dolar_cf <- d$dolar-mean(d$dolar[smf])
d$DxVf <- d$D2011*d$VolTC_cf; d$DxDf <- d$D2011*d$dolar_cf
rob <- rbind(rob, rob_row("Fin de trimestre", rp_fdt ~ VolTC_cf+dolar_cf+D2011+DxVf+DxDf, d, c("DxVf","DxDf"), "D2011"))
# ponderada (tiene hueco 2006T4; na.omit lo colapsa == newey_nogap)
smp <- with(d, !is.na(rp_pond) & !is.na(VolTC) & !is.na(dolar))
d$VolTC_cp <- d$VolTC-mean(d$VolTC[smp]); d$dolar_cp <- d$dolar-mean(d$dolar[smp])
d$DxVp <- d$D2011*d$VolTC_cp; d$DxDp <- d$D2011*d$dolar_cp
rob <- rbind(rob, rob_row("Serie ponderada", rp_pond ~ VolTC_cp+dolar_cp+D2011+DxVp+DxDp, d, c("DxVp","DxDp"), "D2011"))
RES$robustez <- rob
cat("\n== Robustez (nivel y Wald de pendientes) ==\n"); print(rob, digits=4, row.names=FALSE)

## ---- 11. Figuras (ggplot2) -------------------------------------------------
d$fecha_q <- d$anio + (d$trimestre-1)/4
xbreak <- 2011 + (2-1)/4
ggsave("output_R/figuras/Figura_A1_DR.png",
  ggplot(d, aes(fecha_q, DR)) + geom_line(color="navy", linewidth=0.8) +
    geom_hline(yintercept=0, linetype="dashed", color="grey40") +
    geom_vline(xintercept=xbreak, color="red", linewidth=0.8) +
    labs(title="Diferencial de retorno trimestral DR (Gs vs. USD)",
         subtitle="2004-2024; linea roja = Metas de Inflacion (2011:T2)",
         x=NULL, y="DR (pp, trimestral)") + theme_minimal(),
  width=9, height=5, dpi=150)
fis_df <- rbind(
  data.frame(periodo="Pre-2011:T2", comp=c("Prima inflacion","Prima riesgo"),
             valor=c(mean(d$InflDiff[d$D2011==0 & !is.na(d$rp)]), mean(d$rp[d$D2011==0 & !is.na(d$rp)]))),
  data.frame(periodo="Post-2011:T2", comp=c("Prima inflacion","Prima riesgo"),
             valor=c(mean(d$InflDiff[d$D2011==1]), mean(d$rp[d$D2011==1]))))
ggsave("output_R/figuras/Figura_A2_Fisher.png",
  ggplot(fis_df, aes(periodo, valor, fill=comp)) +
    geom_col(position="dodge") + geom_hline(yintercept=0, color="grey40") +
    labs(title="Descomposicion de Fisher del diferencial nominal",
         x=NULL, y="pp anuales", fill=NULL) + theme_minimal(),
  width=8, height=5, dpi=150)
cat("\n[OK] Figuras en output_R/figuras/\n")

## ---- 12. CONTROL AUTOMATICO R vs Stata/Word --------------------------------
esperado <- data.frame(
  concepto = c("Nivel D2011 (principal)","Wald pendientes p","Chow F",
               "Delta rp (Fisher)","Delta InflDiff","Fama beta total","Fama p(beta=1)",
               "Nivel 2011:T1","Wald p 2011:T1","Nivel fin de trimestre","Wald p fin de trim."),
  valor_esperado = c(3.887, 0.142, 26.00, 3.894, -3.065, 1.163, 0.898,
                     3.714, 0.717, 4.023, 0.477))
obtenido <- c(coef(m0$m)["D2011"], w_slopes["p"], w_chow["F"],
              fis$delta[fis$componente=="rp"], fis$delta[fis$componente=="InflDiff"],
              fama["Total","beta"], fama["Total","p_b1"],
              rob$nivel[rob$especificacion=="Sensib. 2011:T1"], rob$Wald_p[rob$especificacion=="Sensib. 2011:T1"],
              rob$nivel[rob$especificacion=="Fin de trimestre"], rob$Wald_p[rob$especificacion=="Fin de trimestre"])
control <- data.frame(esperado, obtenido = as.numeric(obtenido),
                      dif = abs(as.numeric(obtenido) - esperado$valor_esperado))
control$estado <- ifelse(control$dif <= 0.02, "OK", "REVISAR")
cat("\n===================== CONTROL R vs STATA/WORD =====================\n")
print(control, digits=4, row.names=FALSE)
cat(sprintf("\n>> %d/%d cifras coinciden con la tesis (tol 0.02).\n",
            sum(control$estado=="OK"), nrow(control)))

## ---- 13. Exportar todo -----------------------------------------------------
writexl::write_xlsx(c(RES, list(control_R_vs_Stata = control)),
                    "output_R/resultados_tesis_R.xlsx")
cat("\n[OK] Resultados exportados a output_R/resultados_tesis_R.xlsx\n")
cat("\nVEREDICTO (honesto, sin ajustes):\n",
    "- Capa 1 (Fisher): la prima de inflacion cae (~-3.07) y emerge una prima de riesgo\n",
    "  positiva (salto ~+3.89, p<0.001). CONFIRMADA.\n",
    "- Capa 2 (recomposicion interna): la prueba conjunta de pendientes NO rechaza\n",
    "  estabilidad (F=2.00, p=0.142). NO CONCLUYENTE. Se reporta tal cual.\n")
# ============================================================================
