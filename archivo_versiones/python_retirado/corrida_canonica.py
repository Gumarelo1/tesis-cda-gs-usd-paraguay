# -*- coding: utf-8 -*-
"""
CORRIDA CANONICA UNICA - Tesis CDA Gs/USD (Paraguay 2004-2024)
===============================================================
Produce TODOS los numeros del documento final desde UNA sola ejecucion,
con la especificacion CONGELADA tras la auditoria integral del tutor (17/07/2026):

  * Serie primaria: tramo CDA <=365 dias (i_Gs, i_USD), fuente unica SB/BCP.
  * Nucleo: rp = (i_Gs - i_USD) - InflDiff  ("diferencial real residual",
    proxy de prima de riesgo; NO prima pura - lenguaje M4 del informe).
  * D2011 = 1 desde 2011:T1 (PRINCIPAL, especificacion canonica del documento);
    sensibilidad = T2, trimestre del acto institucional (Res.22 Acta 31,
    18-may-2011) [hallazgo M3].
  * Modelo nucleo (Tabla central): rp ~ VolTC_c + Iliq_c + D2011 + DxV + DxI,
    HAC Newey-West lag 4, variables centradas en la muestra del modelo.
  * Sensibilidades: T1, DRcons, L.Iliq, regimen gradual, crisis 08-09,
    HAC lags 3/5/6, post-2011, encaje post-2011, serie ponderada, fdt,
    conversion compuesta (1+i)^(1/4)-1 (hallazgo M7).

La inferencia HAC replica EXACTAMENTE a Stata `newey`: kernel Bartlett,
ajuste de muestra chica T/(T-k), distribucion t(N-k); los Wald usan
F = (Rb)'(R V R')^-1 (Rb)/q ~ F(q, N-k) como `test` de Stata.

Salidas: output/canonica/  (tablas CSV + resultados.json + log.txt + figuras)
"""
import json
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(BASE_DIR)
OUT = os.path.join("output", "canonica")
FIG = os.path.join(OUT, "figuras")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

HAC = 4
R = {}          # JSON maestro de resultados
_loglines = []

def log(*args):
    s = " ".join(str(a) for a in args)
    print(s)
    _loglines.append(s)

# ============================================================================
# newey estilo Stata (kernel Bartlett + ajuste T/(T-k) + t(N-k))
# ============================================================================
def stata_newey(y, X, lag=HAC):
    """y: 1d array; X: 2d (sin constante; se agrega). Devuelve dict."""
    y = np.asarray(y, dtype=float)
    X = np.column_stack([np.asarray(X, dtype=float), np.ones(len(y))]) if X is not None and np.size(X) else np.ones((len(y), 1))
    n, k = X.shape
    XtX_inv = np.linalg.inv(X.T @ X)
    b = XtX_inv @ X.T @ y
    e = y - X @ b
    # HAC Bartlett
    S = (X * e[:, None]).T @ (X * e[:, None])
    for l in range(1, lag + 1):
        w = 1.0 - l / (lag + 1.0)
        G = (X[l:] * e[l:, None]).T @ (X[:-l] * e[:-l, None])
        S += w * (G + G.T)
    S *= n / (n - k)                      # ajuste de muestra chica de Stata
    V = XtX_inv @ S @ XtX_inv
    se = np.sqrt(np.diag(V))
    t = b / se
    df = n - k
    p = 2 * stats.t.sf(np.abs(t), df)
    ss_res = float(e @ e)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return {"b": b, "se": se, "t": t, "p": p, "V": V, "n": n, "k": k,
            "df": df, "r2": r2, "resid": e, "X": X}

def wald(res, idx):
    """Wald F para H0: b[idx]=0 (conjunto), estilo Stata test tras newey."""
    idx = list(idx)
    q = len(idx)
    Rb = res["b"][idx]
    RVR = res["V"][np.ix_(idx, idx)]
    F = float(Rb @ np.linalg.inv(RVR) @ Rb) / q
    p = stats.f.sf(F, q, res["df"])
    return F, p

def lincom(res, w):
    """combinacion lineal w'b con SE HAC -> (est, se, t, p)."""
    w = np.asarray(w, dtype=float)
    est = float(w @ res["b"])
    se = float(np.sqrt(w @ res["V"] @ w))
    t = est / se
    p = 2 * stats.t.sf(abs(t), res["df"])
    return est, se, t, p

def design(df, ycol, xcols, mask=None):
    d = df.dropna(subset=[ycol] + list(xcols)) if xcols else df.dropna(subset=[ycol])
    if mask is not None:
        d = d[mask.loc[d.index]]
    y = d[ycol].values
    X = d[list(xcols)].values if xcols else None
    return y, X, d

# ============================================================================
# 0. CARGA + ASSERTS DE INTEGRIDAD
# ============================================================================
log("=" * 78)
log(" CORRIDA CANONICA - especificacion congelada (auditoria 17/07/2026)")
log("=" * 78)
df = pd.read_csv(os.path.join("data", "processed", "base_cda.csv"))
assert df.shape[0] == 84, f"base debe tener 84 filas, tiene {df.shape[0]}"
per = (df["anio"] * 10 + df["trimestre"]).tolist()
assert per == [a * 10 + q for a in range(2004, 2025) for q in range(1, 5)], "trimestres no continuos"
assert df["i_Gs"].isna().sum() == 0 and df["i_USD"].isna().sum() == 0, "serie <=365 debe estar completa"
assert df["i_Gs_pond"].isna().sum() == 1, "ponderada debe tener 1 faltante (2006Q4)"
assert df["TC"].isna().sum() == 0 and df["dolariz"].isna().sum() == 0
log(f"[OK] Integridad: 84 trimestres 2004Q1-2024Q4; <=365 completa; ponderada 1 faltante; TC/dolariz completas")

# ============================================================================
# 1. VARIABLES
# ============================================================================
df["difn"] = df["i_Gs"] - df["i_USD"]
df["difq"] = df["difn"] / 4.0
df["difq_comp"] = ((1 + df["i_Gs"]/100) ** 0.25 - (1 + df["i_USD"]/100) ** 0.25) * 100  # M7
df["dep"] = df["TC"].pct_change() * 100
df["depLead"] = df["TC"].shift(-1) / df["TC"] * 100 - 100
df["DR"] = df["difq"] - df["depLead"]
df["DRcons"] = df["difq"] - df["depLead"].clip(lower=0)
df["DR_comp"] = df["difq_comp"] - df["depLead"]
df["infl_py"] = (df["IPC_PY"] / df["IPC_PY"].shift(4) - 1) * 100
df["infl_us"] = (df["IPC_US"] / df["IPC_US"].shift(4) - 1) * 100
df["InflDiff"] = df["infl_py"] - df["infl_us"]
df["VolTC"] = df["dep"].rolling(8, min_periods=2).std(ddof=1)
df["Iliq"] = df["dolariz"]
df["rp"] = df["difn"] - df["InflDiff"]
# fdt y ponderada
df["difn_fdt"] = df["i_Gs_fdt"] - df["i_USD_fdt"]
df["difq_fdt"] = df["difn_fdt"] / 4.0
df["difn_pond"] = df["i_Gs_pond"] - df["i_USD_pond"]
df["rp_pond"] = df["difn_pond"] - df["InflDiff"]

# --- regimen: PRINCIPAL desde 2011T1 (especificacion canonica del do-file y
#     del documento: anio calendario del cambio de regimen; justificada en el
#     texto). Sensibilidad: T2 (2011:T2 = trimestre del acto institucional,
#     Res.22 Acta 31 del 18-may-2011) [hallazgo M3 de la auditoria].
df["D2011"]   = (df["anio"] >= 2011).astype(int)
df["D2011t2"] = ((df["anio"] > 2011) | ((df["anio"] == 2011) & (df["trimestre"] >= 2))).astype(int)
log(f"[SPEC] D2011 PRINCIPAL desde 2011T1 (pre n={int((1-df['D2011']).sum())}, post n={int(df['D2011'].sum())}); sensibilidad T2")

MODEL_VARS = ["rp", "VolTC", "Iliq"]
msk_model = df[MODEL_VARS].notna().all(axis=1)
log(f"[SPEC] Muestra del modelo: {msk_model.sum()} trimestres ({df.loc[msk_model,'fecha'].iloc[0]}..{df.loc[msk_model,'fecha'].iloc[-1]})")

# centrado EN LA MUESTRA DEL MODELO
mV = df.loc[msk_model, "VolTC"].mean()
mI = df.loc[msk_model, "Iliq"].mean()
df["VolTC_c"] = df["VolTC"] - mV
df["Iliq_c"] = df["Iliq"] - mI
df["DxV"] = df["D2011"] * df["VolTC_c"]
df["DxI"] = df["D2011"] * df["Iliq_c"]
R["centrado"] = {"media_VolTC": mV, "media_Iliq": mI}

def run_core(dvar, dep="rp", label=""):
    """modelo nucleo con dummy dvar; devuelve res + wald + guarda en R."""
    d = df.copy()
    d["_D"] = d[dvar]
    d["_DxV"] = d["_D"] * d["VolTC_c"]
    d["_DxI"] = d["_D"] * d["Iliq_c"]
    cols = ["VolTC_c", "Iliq_c", "_D", "_DxV", "_DxI"]
    y, X, dd = design(d, dep, cols)
    res = stata_newey(y, X)
    Fp, pp_ = wald(res, [3, 4])         # pendientes: _DxV, _DxI
    Fc, pc = wald(res, [2, 3, 4])       # Chow: _D + interacciones
    out = {"n": res["n"], "r2": res["r2"],
           "b": dict(zip(["VolTC_c","Iliq_c","D2011","DxV","DxI","cons"], res["b"].tolist())),
           "se": dict(zip(["VolTC_c","Iliq_c","D2011","DxV","DxI","cons"], res["se"].tolist())),
           "p": dict(zip(["VolTC_c","Iliq_c","D2011","DxV","DxI","cons"], res["p"].tolist())),
           "wald_pendientes": {"F": Fp, "p": pp_}, "chow": {"F": Fc, "p": pc}}
    if label:
        R[label] = out
    return res, out

# ============================================================================
# ETAPA 1 - Descriptivos DR / DRcons / rp (cortes con T2)
# ============================================================================
log("\n" + "=" * 78); log(" ETAPA 1 - Descriptivos (corte 2011T2)"); log("=" * 78)
rows = []
for v in ["DR", "DRcons", "rp"]:
    for lbl, m in [("Total", df.index >= 0), ("Pre", df["D2011"] == 0), ("Post", df["D2011"] == 1)]:
        s = df.loc[m, v].dropna()
        rows.append({"serie": v, "corte": lbl, "n": len(s), "media": s.mean(),
                     "mediana": s.median(), "sd": s.std(ddof=1), "min": s.min(), "max": s.max(),
                     "asimetria": float(stats.skew(s, bias=True)),
                     "curtosis": float(stats.kurtosis(s, fisher=False, bias=True)),
                     "pct_pos": 100 * (s > 0).mean()})
desc = pd.DataFrame(rows)
desc.to_csv(os.path.join(OUT, "t1_descriptivos.csv"), index=False)
R["descriptivos"] = desc.to_dict("records")
log(desc.round(3).to_string(index=False))

# ============================================================================
# ETAPA 2 - Fama (difq; sensibilidad compuesta y fdt) + DRcons unilateral
# ============================================================================
log("\n" + "=" * 78); log(" ETAPA 2 - Fama (regresor trimestralizado difq) + DRcons"); log("=" * 78)
fama_rows = []
for reglab, regcol in [("difq", "difq"), ("difq_comp (M7)", "difq_comp"), ("difq_fdt", "difq_fdt")]:
    for lbl, m in [("Total", None), ("Pre", df["D2011"] == 0), ("Post", df["D2011"] == 1)]:
        y, X, dd = design(df, "depLead", [regcol], m)
        res = stata_newey(y, X)
        b, se = res["b"][0], res["se"][0]
        t1 = (b - 1) / se
        p1 = 2 * stats.t.sf(abs(t1), res["df"])
        fama_rows.append({"regresor": reglab, "muestra": lbl, "n": res["n"], "beta": b,
                          "se": se, "p_b0": res["p"][0], "t_b1": t1, "p_b1": p1, "r2": res["r2"]})
fama = pd.DataFrame(fama_rows)
fama.to_csv(os.path.join(OUT, "t2_fama.csv"), index=False)
R["fama"] = fama.to_dict("records")
log(fama.round(4).to_string(index=False))

log("\n-- DRcons: H0 mu<=0 vs HA mu>0 (HAC) --")
cons_rows = []
for lbl, m in [("Total", None), ("Pre", df["D2011"] == 0), ("Post", df["D2011"] == 1)]:
    y, X, dd = design(df, "DRcons", [], m)
    res = stata_newey(y, None)
    mu, se = res["b"][0], res["se"][0]
    t = mu / se
    p1 = stats.t.sf(t, res["df"])
    tc = stats.t.ppf(0.975, res["df"])
    cons_rows.append({"muestra": lbl, "n": res["n"], "media": mu, "se_hac": se, "t": t,
                      "p_unilateral": p1, "ic_low": mu - tc * se, "ic_high": mu + tc * se,
                      "pct_pos": 100 * (dd["DRcons"] > 0).mean()})
cons = pd.DataFrame(cons_rows)
cons.to_csv(os.path.join(OUT, "t2_drcons.csv"), index=False)
R["drcons"] = cons.to_dict("records")
log(cons.round(4).to_string(index=False))

# ============================================================================
# ETAPA 3 - Integracion (ADF+KPSS) + quiebre endogeno + dif. de medias
# ============================================================================
log("\n" + "=" * 78); log(" ETAPA 3 - Integracion + quiebre estructural"); log("=" * 78)
from statsmodels.tsa.stattools import adfuller, kpss
import warnings
integ_rows = []
for v in ["DR", "rp", "VolTC", "Iliq", "InflDiff"]:
    s = df[v].dropna().values
    adf = adfuller(s, autolag="BIC")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        kp = kpss(s, regression="c", nlags=4)
    integ_rows.append({"serie": v, "n": len(s), "ADF": adf[0], "p_ADF": adf[1], "lags_ADF": adf[2],
                       "KPSS": kp[0], "p_KPSS": kp[1],
                       "veredicto": ("I(0)" if (adf[1] < 0.05 and kp[1] > 0.05)
                                     else "I(1)" if (adf[1] >= 0.05 and kp[1] <= 0.05) else "mixto")})
integ = pd.DataFrame(integ_rows)
integ.to_csv(os.path.join(OUT, "t3_integracion.csv"), index=False)
R["integracion"] = integ.to_dict("records")
log(integ.round(4).to_string(index=False))

# quiebre endogeno (Bai-Perron via ruptures, sobre rp y DR)
try:
    import ruptures as rpt
    for v in ["rp", "DR"]:
        d = df.dropna(subset=[v])
        sig = d[v].values.reshape(-1, 1)
        algo = rpt.Binseg(model="l2").fit(sig)
        for nb in [1, 2, 3]:
            bk = algo.predict(n_bkps=nb)[:-1]
            fechas = [d["fecha"].iloc[i] for i in bk]
            if nb == 1:
                R[f"quiebre_{v}"] = fechas[0]
            log(f"  Bai-Perron {v}: {nb} quiebre(s) -> {fechas}")
except Exception as e:
    log("  [aviso] ruptures no disponible:", e)

log("\n-- Dif. de medias pre/post (HAC): rp, DR --")
for v in ["rp", "DR"]:
    y, X, dd = design(df, v, ["D2011"])
    res = stata_newey(y, X)
    R[f"difmedias_{v}"] = {"coef": res["b"][0], "se": res["se"][0], "p": res["p"][0], "n": res["n"]}
    log(f"  {v}: D2011={res['b'][0]:+.3f}  se={res['se'][0]:.3f}  p={res['p'][0]:.4f}  n={res['n']}")

# ============================================================================
# ETAPA 4 - Descomposicion de Fisher (niveles pre/post)
# ============================================================================
log("\n" + "=" * 78); log(" ETAPA 4 - Fisher: difn = InflDiff + rp (niveles, corte T2)"); log("=" * 78)
d4 = df.dropna(subset=["difn", "InflDiff", "rp"])
fis_rows = []
for lbl, m in [("Pre", d4["D2011"] == 0), ("Post", d4["D2011"] == 1)]:
    sub = d4[m]
    fis_rows.append({"corte": lbl, "n": len(sub), "difn": sub["difn"].mean(),
                     "InflDiff": sub["InflDiff"].mean(), "rp": sub["rp"].mean(),
                     "sd_InflDiff": sub["InflDiff"].std(ddof=1), "sd_rp": sub["rp"].std(ddof=1)})
fis = pd.DataFrame(fis_rows)
fis.to_csv(os.path.join(OUT, "t4_fisher.csv"), index=False)
R["fisher"] = fis.to_dict("records")
log(fis.round(3).to_string(index=False))
log("\n-- t HAC de la diferencia pre/post (muestra de la descomposicion, n=80) --")
for v in ["rp", "InflDiff", "difn"]:
    y, X, dd = design(df.dropna(subset=["rp"]), v, ["D2011"])
    res = stata_newey(y, X)
    R[f"fisher_t_{v}"] = {"delta": res["b"][0], "se": res["se"][0], "p": res["p"][0]}
    log(f"  {v}: delta={res['b'][0]:+.3f}  p={res['p'][0]:.4f}")

# ============================================================================
# ETAPA 5 (NUCLEO) - rp centrado, D2011 desde T2
# ============================================================================
log("\n" + "=" * 78); log(" ETAPA 5 (NUCLEO) - rp ~ VolTC_c + Iliq_c + D2011(T2) + DxV + DxI"); log("=" * 78)
res_core, out_core = run_core("D2011", "rp", "modelo_nucleo")
names = ["VolTC_c", "Iliq_c", "D2011", "DxV", "DxI", "constante"]
tbl = pd.DataFrame({"variable": names, "coef": res_core["b"], "se_hac": res_core["se"],
                    "t": res_core["t"], "p": res_core["p"]})
tbl.to_csv(os.path.join(OUT, "t5_modelo_nucleo.csv"), index=False)
log(tbl.round(4).to_string(index=False))
log(f"  n={out_core['n']}  R2={out_core['r2']:.3f}")
log(f"  Wald pendientes (DxV=DxI=0): F={out_core['wald_pendientes']['F']:.3f}  p={out_core['wald_pendientes']['p']:.4f}")
log(f"  Chow (D2011=DxV=DxI=0):     F={out_core['chow']['F']:.3f}  p={out_core['chow']['p']:.4f}")
# unilaterales b5/b6
t5 = res_core["b"][3] / res_core["se"][3]
t6 = res_core["b"][4] / res_core["se"][4]
R["unilaterales"] = {"p_b5_menor0": float(stats.t.cdf(t5, res_core["df"])),
                     "p_b6_mayor0": float(stats.t.sf(t6, res_core["df"]))}
log(f"  b5 (DxV) p_1cola(<0)={R['unilaterales']['p_b5_menor0']:.3f} | b6 (DxI) p_1cola(>0)={R['unilaterales']['p_b6_mayor0']:.3f}")

# dominancia estandarizada
dm = df[msk_model].copy()
dm["zV"] = (dm["VolTC"] - dm["VolTC"].mean()) / dm["VolTC"].std(ddof=1)
dm["zI"] = (dm["Iliq"] - dm["Iliq"].mean()) / dm["Iliq"].std(ddof=1)
dm["DzV"] = dm["D2011"] * dm["zV"]
dm["DzI"] = dm["D2011"] * dm["zI"]
y, X, dd = design(dm, "rp", ["zV", "zI", "D2011", "DzV", "DzI"])
res_z = stata_newey(y, X)
inc_V = lincom(res_z, [1, 0, 0, 1, 0, 0])
inc_I = lincom(res_z, [0, 1, 0, 0, 1, 0])
domin = lincom(res_z, [-1, 1, 0, -1, 1, 0])
R["dominancia"] = {"incidencia_post_VolTC_sd": inc_V[0], "incidencia_post_Iliq_sd": inc_I[0],
                   "dif": domin[0], "t": domin[2], "p": domin[3]}
log(f"\n-- Dominancia post (estandarizada): VolTC {inc_V[0]:+.3f} sd | Iliq {inc_I[0]:+.3f} sd | dif p={domin[3]:.3f}")

# JB de residuos
e = res_core["resid"]
sk = stats.skew(e, bias=True); ku = stats.kurtosis(e, fisher=False, bias=True)
JB = len(e) / 6 * (sk ** 2 + (ku - 3) ** 2 / 4)
R["jarque_bera"] = {"JB": float(JB), "p": float(stats.chi2.sf(JB, 2))}
log(f"-- JB residuos: {JB:.2f} (p={R['jarque_bera']['p']:.3f})")

# Engle-Granger ACOTADO: relacion de largo plazo SIN dummies, CV MacKinnon
from statsmodels.tsa.stattools import coint
d5 = df.dropna(subset=["rp", "VolTC", "Iliq"])
tEG, pEG, cvEG = coint(d5["rp"], d5[["VolTC", "Iliq"]], trend="c", autolag="BIC")
R["engle_granger"] = {"t": float(tEG), "p": float(pEG), "cv_5pct": float(cvEG[1]),
                      "nota": "relacion de largo plazo SIN dummies; CV de MacKinnon para cointegracion (no ADF ordinarios)"}
log(f"-- Engle-Granger (rp ~ VolTC + Iliq, sin dummies, CV MacKinnon): t={tEG:.3f}  p={pEG:.4f}  cv5%={cvEG[1]:.3f}")

# ============================================================================
# ETAPA 6 - Diagnosticos del OLS paralelo
# ============================================================================
log("\n" + "=" * 78); log(" ETAPA 6 - Diagnosticos"); log("=" * 78)
import statsmodels.api as sm
from statsmodels.stats.diagnostic import acorr_breusch_godfrey, het_white, linear_reset
from statsmodels.stats.outliers_influence import variance_inflation_factor
dd6 = df.dropna(subset=["rp", "VolTC_c", "Iliq_c"]).copy()
X6 = sm.add_constant(dd6[["VolTC_c", "Iliq_c", "D2011", "DxV", "DxI"]])
ols6 = sm.OLS(dd6["rp"], X6).fit()
diag = {}
for L in [1, 2, 3, 4]:
    bg = acorr_breusch_godfrey(ols6, nlags=L)
    diag[f"BG_lag{L}"] = {"LM": bg[0], "p": bg[1]}
w = het_white(ols6.resid, X6)
diag["White"] = {"LM": w[0], "p": w[1]}
rst = linear_reset(ols6, power=3, use_f=True)
diag["RESET"] = {"F": float(rst.fvalue), "p": float(rst.pvalue)}
vifs = {X6.columns[i]: float(variance_inflation_factor(X6.values, i)) for i in range(1, X6.shape[1])}
diag["VIF"] = vifs
R["diagnosticos"] = diag
for k, v in diag.items():
    log(f"  {k}: {v}")

# ============================================================================
# ROBUSTEZ (todas en un solo bloque, tabla resumen)
# ============================================================================
log("\n" + "=" * 78); log(" ROBUSTEZ"); log("=" * 78)
rob_rows = []

def add_rob(nombre, out, extra=""):
    rob_rows.append({"especificacion": nombre, "n": out["n"],
                     "D2011": out["b"]["D2011"], "p_D2011": out["p"]["D2011"],
                     "DxV": out["b"]["DxV"], "p_DxV": out["p"]["DxV"],
                     "DxI": out["b"]["DxI"], "p_DxI": out["p"]["DxI"],
                     "Wald_F": out["wald_pendientes"]["F"], "Wald_p": out["wald_pendientes"]["p"],
                     "nota": extra})

add_rob("PRINCIPAL: rp, D2011 desde T1", out_core)

# (a) sensibilidad T2 (trimestre del acto institucional; M3)
_, out_t2 = run_core("D2011t2", "rp", "modelo_T2")
add_rob("Sensibilidad: D2011 desde T2 (M3)", out_t2)

# (b) DRcons como dependiente
_, out_dc = run_core("D2011", "DRcons", "modelo_DRcons")
add_rob("DRcons como dependiente", out_dc)

# (c) dolarizacion rezagada
dl = df.copy()
dl["LIliq"] = dl["Iliq"].shift(1)
mM = dl[["rp", "VolTC", "LIliq"]].notna().all(axis=1)
dl["LIliq_c"] = dl["LIliq"] - dl.loc[mM, "LIliq"].mean()
dl["DxLI"] = dl["D2011"] * dl["LIliq_c"]
y, X, _ = design(dl, "rp", ["VolTC_c", "LIliq_c", "D2011", "DxV", "DxLI"])
res_l = stata_newey(y, X)
Fl, pl = wald(res_l, [3, 4])
rob_rows.append({"especificacion": "Iliq rezagada (t-1)", "n": res_l["n"],
                 "D2011": res_l["b"][2], "p_D2011": res_l["p"][2],
                 "DxV": res_l["b"][3], "p_DxV": res_l["p"][3],
                 "DxI": res_l["b"][4], "p_DxI": res_l["p"][4],
                 "Wald_F": Fl, "Wald_p": pl, "nota": "simultaneidad"})

# (d) regimen gradual (rampa 2011T2-2016, 1 desde 2017)
dg = df.copy()
idx_ramp = (dg["anio"] * 4 + dg["trimestre"]) - (2011 * 4 + 1) + 1   # rampa desde 2011T1 (como el do-file)
dg["Reg"] = np.clip(idx_ramp / 24.0, 0, None)
dg.loc[dg["anio"] < 2011, "Reg"] = 0
dg["Reg"] = dg["Reg"].clip(upper=1)
dg["RxV"] = dg["Reg"] * dg["VolTC_c"]
dg["RxI"] = dg["Reg"] * dg["Iliq_c"]
y, X, _ = design(dg, "rp", ["VolTC_c", "Iliq_c", "Reg", "RxV", "RxI"])
res_g = stata_newey(y, X)
Fg, pg = wald(res_g, [3, 4])
rob_rows.append({"especificacion": "Regimen gradual (rampa)", "n": res_g["n"],
                 "D2011": res_g["b"][2], "p_D2011": res_g["p"][2],
                 "DxV": res_g["b"][3], "p_DxV": res_g["p"][3],
                 "DxI": res_g["b"][4], "p_DxI": res_g["p"][4],
                 "Wald_F": Fg, "Wald_p": pg, "nota": "intensidad 0-1"})

# (e) control crisis 2008-09
dc = df.copy()
dc["Crisis"] = (((dc["anio"] == 2008) & (dc["trimestre"] >= 2)) | ((dc["anio"] == 2009) & (dc["trimestre"] <= 3))).astype(int)
y, X, _ = design(dc, "rp", ["VolTC_c", "Iliq_c", "D2011", "DxV", "DxI", "Crisis"])
res_c = stata_newey(y, X)
Fc2, pc2 = wald(res_c, [3, 4])
rob_rows.append({"especificacion": "Control crisis 2008Q2-2009Q3", "n": res_c["n"],
                 "D2011": res_c["b"][2], "p_D2011": res_c["p"][2],
                 "DxV": res_c["b"][3], "p_DxV": res_c["p"][3],
                 "DxI": res_c["b"][4], "p_DxI": res_c["p"][4],
                 "Wald_F": Fc2, "Wald_p": pc2,
                 "nota": f"Crisis={res_c['b'][5]:+.2f} (p={res_c['p'][5]:.3f})"})
R["crisis_control"] = {"D2011": res_c["b"][2], "p": res_c["p"][2],
                       "crisis_coef": res_c["b"][5], "crisis_p": res_c["p"][5]}

# (f) HAC lags
for L in [3, 5, 6]:
    y, X, _ = design(df, "rp", ["VolTC_c", "Iliq_c", "D2011", "DxV", "DxI"])
    res_h = stata_newey(y, X, lag=L)
    Fh, ph = wald(res_h, [3, 4])
    rob_rows.append({"especificacion": f"HAC lag={L}", "n": res_h["n"],
                     "D2011": res_h["b"][2], "p_D2011": res_h["p"][2],
                     "DxV": res_h["b"][3], "p_DxV": res_h["p"][3],
                     "DxI": res_h["b"][4], "p_DxI": res_h["p"][4],
                     "Wald_F": Fh, "Wald_p": ph, "nota": ""})

# (g) serie ponderada (sensibilidad de definicion)
y, X, _ = design(df.assign(dep_p=df["rp_pond"]), "dep_p", ["VolTC_c", "Iliq_c", "D2011", "DxV", "DxI"])
res_p = stata_newey(y, X)
Fp2, pp2 = wald(res_p, [3, 4])
rob_rows.append({"especificacion": "Serie PONDERADA (sensibilidad)", "n": res_p["n"],
                 "D2011": res_p["b"][2], "p_D2011": res_p["p"][2],
                 "DxV": res_p["b"][3], "p_DxV": res_p["p"][3],
                 "DxI": res_p["b"][4], "p_DxI": res_p["p"][4],
                 "Wald_F": Fp2, "Wald_p": pp2, "nota": "empalme 2010/11 + composicion de plazos"})
R["ponderada"] = {"D2011": res_p["b"][2], "p_D2011": res_p["p"][2], "Wald_F": Fp2, "Wald_p": pp2, "n": res_p["n"]}

# (h) post-2011 submuestra (sin dummies)
y, X, _ = design(df, "rp", ["VolTC", "Iliq"], df["D2011"] == 1)
res_post = stata_newey(y, X)
R["post2011"] = {"n": res_post["n"], "VolTC": res_post["b"][0], "p_VolTC": res_post["p"][0],
                 "Iliq": res_post["b"][1], "p_Iliq": res_post["p"][1], "r2": res_post["r2"]}
log(f"\n-- Post-2011 (sin dummies): VolTC={res_post['b'][0]:+.3f} (p={res_post['p'][0]:.3f}), Iliq={res_post['b'][1]:+.3f} (p={res_post['p'][1]:.3f}), n={res_post['n']}")

# (i) encaje post-2011 (efectivo)
if df["encaje_me_efec"].notna().sum() > 10:
    de = df.copy()
    de["brecha"] = de["encaje_me_efec"] - de["encaje_mn_efec"]
    sub = de[(de["D2011"] == 1)].dropna(subset=["rp", "VolTC", "brecha"])
    if sub["brecha"].std(ddof=1) > 0:
        y, X, _ = design(sub, "rp", ["VolTC", "brecha"])
        res_e = stata_newey(y, X)
        R["encaje_post"] = {"n": res_e["n"], "brecha": res_e["b"][1], "p": res_e["p"][1]}
        log(f"-- Encaje post-2011 (brecha efectiva): coef={res_e['b'][1]:+.3f} (p={res_e['p'][1]:.3f}), n={res_e['n']}")

rob = pd.DataFrame(rob_rows)
rob.to_csv(os.path.join(OUT, "t7_robustez.csv"), index=False)
log("\n" + rob.round(3).to_string(index=False))

# ============================================================================
# COMPOSICION INTERNA (Shapley) + ventana movil (figura)
# ============================================================================
log("\n" + "=" * 78); log(" Composicion interna de rp (Shapley) + ventana movil"); log("=" * 78)
shap_rows = []
for lbl, m in [("Pre", df["D2011"] == 0), ("Post", df["D2011"] == 1)]:
    sub = df[m].dropna(subset=["rp", "VolTC", "Iliq"])
    r2v = sm.OLS(sub["rp"], sm.add_constant(sub[["VolTC"]])).fit().rsquared
    r2i = sm.OLS(sub["rp"], sm.add_constant(sub[["Iliq"]])).fit().rsquared
    r2b = sm.OLS(sub["rp"], sm.add_constant(sub[["VolTC", "Iliq"]])).fit().rsquared
    cv = 0.5 * (r2v + r2b - r2i); ci = 0.5 * (r2i + r2b - r2v)
    shap_rows.append({"corte": lbl, "n": len(sub), "R2": r2b,
                      "pct_VolTC": 100 * cv / (cv + ci), "pct_Iliq": 100 * ci / (cv + ci)})
shap = pd.DataFrame(shap_rows)
shap.to_csv(os.path.join(OUT, "t8_shapley.csv"), index=False)
R["shapley"] = shap.to_dict("records")
log(shap.round(3).to_string(index=False))

# ventana movil 20T
dwin = df.dropna(subset=["rp", "VolTC", "Iliq"]).reset_index(drop=True)
W = 20
mov = []
for i in range(W - 1, len(dwin)):
    sub = dwin.iloc[i - W + 1:i + 1]
    r2v = sm.OLS(sub["rp"], sm.add_constant(sub[["VolTC"]])).fit().rsquared
    r2i = sm.OLS(sub["rp"], sm.add_constant(sub[["Iliq"]])).fit().rsquared
    r2b = sm.OLS(sub["rp"], sm.add_constant(sub[["VolTC", "Iliq"]])).fit().rsquared
    cv = 0.5 * (r2v + r2b - r2i); ci = 0.5 * (r2i + r2b - r2v)
    mov.append({"fecha": sub["fecha"].iloc[-1], "pct_VolTC": 100 * cv / (cv + ci),
                "pct_Iliq": 100 * ci / (cv + ci), "R2": r2b})
movdf = pd.DataFrame(mov)
movdf.to_csv(os.path.join(OUT, "t9_ventana_movil.csv"), index=False)

# ============================================================================
# FIGURAS
# ============================================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
tt = np.arange(len(df))
lbls = df["fecha"].values
brk = int(df.index[(df["anio"] == 2011) & (df["trimestre"] == 1)][0])

def xticks(ax):
    ix = list(range(0, len(df), 8))
    ax.set_xticks(ix); ax.set_xticklabels([lbls[i] for i in ix], rotation=45, fontsize=8)

# F1: rp
fig, ax = plt.subplots(figsize=(11, 5.2))
ax.plot(tt, df["rp"], color="#1f4e78", lw=1.8)
ax.axhline(0, ls=":", color="gray"); ax.axvline(brk, color="crimson", lw=1.6)
ax.set_title("Diferencial real residual (rp) — marca: régimen MdI desde 2011")
ax.set_ylabel("pp anuales"); xticks(ax); fig.tight_layout()
fig.savefig(os.path.join(FIG, "f1_rp_serie.png"), dpi=150); plt.close(fig)

# F2: DR y DRcons
fig, ax = plt.subplots(figsize=(11, 5.2))
ax.plot(tt, df["DR"], color="#1f4e78", lw=1.6, label="DR (estándar)")
ax.plot(tt, df["DRcons"], color="#e07b39", lw=1.4, ls="--", label="DRcons (conservador)")
ax.axhline(0, ls=":", color="gray"); ax.axvline(brk, color="crimson", lw=1.6)
ax.legend(); ax.set_title("Diferencial de retorno trimestral DR y DRcons")
ax.set_ylabel("pp trimestrales"); xticks(ax); fig.tight_layout()
fig.savefig(os.path.join(FIG, "f2_DR_DRcons.png"), dpi=150); plt.close(fig)

# F3: Fisher (barras pre/post)
fig, ax = plt.subplots(figsize=(8, 5))
pre = fis.iloc[0]; post = fis.iloc[1]
x = np.arange(2)
ax.bar(x - 0.18, [pre["InflDiff"], post["InflDiff"]], width=0.36, label="Prima de inflación (InflDiff)", color="#7a9cc6")
ax.bar(x + 0.18, [pre["rp"], post["rp"]], width=0.36, label="Diferencial real residual (rp)", color="#1f4e78")
for xi, tot in zip(x, [pre["difn"], post["difn"]]):
    ax.text(xi, max(pre["difn"], post["difn"]) * 1.02, f"difn={tot:.2f}", ha="center", fontsize=9, fontweight="bold")
ax.axhline(0, color="gray", lw=0.8)
ax.set_xticks(x); ax.set_xticklabels([f"Pre-2011 (n={int(pre['n'])})", f"Post (n={int(post['n'])})"])
ax.set_ylabel("pp anuales"); ax.legend(); ax.set_title("Descomposición de Fisher del diferencial nominal")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "f3_fisher.png"), dpi=150); plt.close(fig)

# F4: composicion movil
fig, ax = plt.subplots(figsize=(11, 5.2))
tw = np.arange(len(movdf))
ax.plot(tw, movdf["pct_VolTC"], color="#a4262c", lw=1.8, label="Riesgo cambiario (VolTC)")
ax.plot(tw, movdf["pct_Iliq"], color="#1f4e78", lw=1.8, label="Dolarización/segmentación (Iliq)")
ax.axhline(50, ls=":", color="gray")
bw = movdf.index[movdf["fecha"] == "2011Q1"]
if len(bw): ax.axvline(int(bw[0]), color="black", ls="--")
ix = list(range(0, len(movdf), 8))
ax.set_xticks(ix); ax.set_xticklabels([movdf["fecha"].iloc[i] for i in ix], rotation=45, fontsize=8)
ax.set_ylabel("% del R² explicado"); ax.legend()
ax.set_title("Composición del diferencial real residual (ventana móvil de 20 trimestres)")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "f4_composicion_movil.png"), dpi=150); plt.close(fig)

# F5: TPM vs i_Gs (descriptivo post-2011)
dpm = df.dropna(subset=["tpm"])
fig, ax = plt.subplots(figsize=(11, 5.2))
t5x = np.arange(len(dpm))
ax.plot(t5x, dpm["i_Gs"], color="#1f4e78", lw=1.8, label="CDA ≤365 días en Gs (tasa efectiva)")
ax.plot(t5x, dpm["tpm"], color="#a4262c", lw=1.8, ls="--", label="TPM (BCP)")
ix = list(range(0, len(dpm), 6))
ax.set_xticks(ix); ax.set_xticklabels([dpm["fecha"].iloc[i] for i in ix], rotation=45, fontsize=8)
ax.set_ylabel("% anual"); ax.legend(); ax.set_title("Tasa CDA en guaraníes y TPM (post-2011)")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "f5_tpm_iGs.png"), dpi=150); plt.close(fig)
corr_tpm = dpm["i_Gs"].corr(dpm["tpm"])
spread_tpm = (dpm["i_Gs"] - dpm["tpm"]).mean()
R["tpm"] = {"n": int(len(dpm)), "corr_iGs": float(corr_tpm), "spread_medio": float(spread_tpm)}
log(f"\n-- TPM: corr(i_Gs, TPM) post-2011 = {corr_tpm:.3f}; spread medio i_Gs-TPM = {spread_tpm:+.2f} pp (n={len(dpm)})")

# ============================================================================
# GUARDAR
# ============================================================================
with open(os.path.join(OUT, "resultados.json"), "w", encoding="utf-8") as f:
    json.dump(R, f, indent=1, ensure_ascii=False, default=float)
with open(os.path.join(OUT, "log.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(_loglines))
log("\n[FIN] Salidas en output/canonica/ (tablas, figuras, resultados.json, log.txt)")
