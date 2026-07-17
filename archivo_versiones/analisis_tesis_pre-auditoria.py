# -*- coding: utf-8 -*-
"""
================================================================================
 TESIS: Diferencial de Retornos CDA Guaranies vs Dolares
        Sistema Bancario Paraguayo, 2004-2024  (Scavone & Fernandez)
--------------------------------------------------------------------------------
 SCRIPT PYTHON  ==  MOTOR DE EJECUCION Y VALIDACION CRUZADA
 Equivalente 1:1 al do-file de Stata 17:  stata/tesis_cda.do
--------------------------------------------------------------------------------
 Implementa las 6 etapas del marco metodologico (Seccion 7.2 del anteproyecto).
 Genera todas las tablas (output/tablas/*.csv), figuras (output/figuras/*.png)
 y un JSON de coincidencia para el cross-check contra Stata.

 HONESTIDAD DE RESULTADOS: los coeficientes, SE, p-valores e IC se reportan
 TAL COMO SALEN. No se ajusta ningun dato ni especificacion para forzar la
 hipotesis. Un resultado nulo o contrario, bien documentado, es valido.

 Requisitos: pandas, numpy, statsmodels, scipy, ruptures, matplotlib.
================================================================================
"""

import os
import json
import warnings

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import statsmodels.api as sm
from statsmodels.stats.diagnostic import (
    acorr_breusch_godfrey,
    het_white,
    linear_reset,
)
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tsa.stattools import adfuller
from scipy import stats

try:
    import ruptures as rpt

    HAS_RUPTURES = True
except Exception:  # pragma: no cover
    HAS_RUPTURES = False

warnings.simplefilter("ignore")  # silenciar avisos de statsmodels sobre muestras chicas

# ------------------------------------------------------------------------------
# 0. CONFIGURACION Y PATHS
# ------------------------------------------------------------------------------
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "processed", "base_cda.csv")
OUT = os.path.join(ROOT, "output")
TAB = os.path.join(OUT, "tablas")
FIG = os.path.join(OUT, "figuras")
for d in (OUT, TAB, FIG):
    os.makedirs(d, exist_ok=True)

# Bandwidth HAC (Newey-West).  Andrews(1991) O(T^{1/3}); regla Newey-West
# floor(4*(T/100)^(2/9)).  Con T=84 ambas dan ~4.  Se documenta y se fija LAG=4.
HAC_LAG = 4
YEAR_BREAK = 2011  # Metas de Inflacion (regimen del BCP)

# almacen de numeros clave para el cross-check con Stata
XCHECK = {}


def q_index(anio, trimestre):
    """Indice trimestral entero continuo (equivalente a yq() de Stata)."""
    return (anio - 2004) * 4 + (trimestre - 1)


def fmt(x, nd=4):
    if x is None or (isinstance(x, float) and (np.isnan(x))):
        return "."
    return f"{x:.{nd}f}"


# ==============================================================================
# CARGA Y CONSTRUCCION DE VARIABLES (Etapa 1 - insumos)
# ==============================================================================
def load_and_build():
    """
    Carga base_cda.csv y construye TODAS las variables derivadas.
    Reglas identicas al do-file de Stata (mismas formulas y alineacion temporal).

    Condiciones del auditor respetadas:
      - i_Gs / i_USD vacios en 2006Q4 y 2010Q4 -> NaN (exclusion listwise, sin imputar).
      - TC = promedio mensual del mes de cierre (proxy declarado).
      - Empalme 2010/2011 (~0.7-0.8pp): absorbido por D2011.
    """
    df = pd.read_csv(DATA)

    # indice trimestral continuo + orden
    df["t"] = df.apply(lambda r: q_index(r["anio"], r["trimestre"]), axis=1)
    df = df.sort_values("t").reset_index(drop=True)

    # --- diferencial de tasas ---
    df["difn"] = df["i_Gs"] - df["i_USD"]          # nominal (pp, % anual)
    df["difq"] = df["difn"] / 4.0                   # trimestralizado

    # --- depreciacion del guarani (%) ---
    df["dep"] = (df["TC"] - df["TC"].shift(1)) / df["TC"].shift(1) * 100.0        # dep_t
    df["depLead"] = (df["TC"].shift(-1) - df["TC"]) / df["TC"] * 100.0            # dep_{t+1} (= F.dep)

    # --- series de retorno (alineadas en t con determinantes en t) ---
    df["DR"] = df["difq"] - df["depLead"]                        # estandar (= DR_{t+1} de la tesis)
    df["DRcons"] = df["difq"] - df["depLead"].clip(lower=0)      # conservadora: max(depLead,0)

    # --- regimen Metas de Inflacion ---
    df["D2011"] = (df["anio"] >= YEAR_BREAK).astype(int)

    # --- diferencial de inflacion interanual (Fisher) desde indices IPC ---
    df["infl_py"] = (df["IPC_PY"] - df["IPC_PY"].shift(4)) / df["IPC_PY"].shift(4) * 100.0
    df["infl_us"] = (df["IPC_US"] - df["IPC_US"].shift(4)) / df["IPC_US"].shift(4) * 100.0
    df["InflDiff"] = df["infl_py"] - df["infl_us"]

    # --- volatilidad cambiaria: sd movil 8 trim de la depreciacion ---
    # Replica EXACTA de:  rangestat (sd) VolTC = dep, interval(t -7 0)
    # rangestat calcula sd (ddof=1) sobre las obs NO-missing de la ventana calendario
    # [t-7, t], exigiendo n>=2.  dep solo falta en 2004Q1, por lo que la unica
    # diferencia frente a min_periods=8 esta en el borde inicial (2004Q3-2005Q4).
    df["VolTC"] = df["dep"].rolling(window=8, min_periods=2).std(ddof=1)

    # --- prima de riesgo residual (Etapa 4) ---
    df["rp"] = df["difn"] - df["InflDiff"]

    # --- proxy de iliquidez ---
    df["Iliq"] = df["dolariz"]

    # --- interacciones (Etapa 5) ---
    df["D_VolTC"] = df["D2011"] * df["VolTC"]
    df["D_Iliq"] = df["D2011"] * df["Iliq"]

    # variantes fin-de-trimestre (robustez)
    df["difn_fdt"] = df["i_Gs_fdt"] - df["i_USD_fdt"]
    df["difq_fdt"] = df["difn_fdt"] / 4.0
    df["DR_fdt"] = df["difq_fdt"] - df["depLead"]
    df["DRcons_fdt"] = df["difq_fdt"] - df["depLead"].clip(lower=0)

    return df


# ==============================================================================
# UTILIDADES ECONOMETRICAS
# ==============================================================================
def describe(series, label):
    """Estadistica descriptiva completa (media, mediana, sd, CV, min, max, skew, kurt)."""
    s = series.dropna()
    n = len(s)
    mean = s.mean()
    sd = s.std(ddof=1)
    cv = sd / mean if mean != 0 else np.nan
    return {
        "serie": label,
        "n": n,
        "media": mean,
        "mediana": s.median(),
        "sd": sd,
        "CV": cv,
        "min": s.min(),
        "max": s.max(),
        # skew/kurt con las mismas convenciones que Stata:
        #   Stata skewness = g1 (sesgada, no corregida por muestra)
        #   Stata kurtosis = kurtosis cruda (no exceso): normal = 3
        "asimetria": stats.skew(s, bias=True),
        "curtosis": stats.kurtosis(s, fisher=False, bias=True),
    }


def newey_ols(y, X, lag=HAC_LAG, add_const=True):
    """
    MCO con errores estandar HAC (Newey-West), replicando `newey ... , lag(#)`.
    statsmodels cov_type='HAC' con maxlags=lag y use_correction=True == Stata newey.
    """
    yv = pd.Series(y).astype(float)
    Xv = pd.DataFrame(X).astype(float)
    d = pd.concat([yv, Xv], axis=1).dropna()
    yy = d.iloc[:, 0]
    XX = d.iloc[:, 1:]
    if add_const:
        XX = sm.add_constant(XX, has_constant="add")
    model = sm.OLS(yy, XX)
    # use_t=True  =>  p-valores e IC con distribucion t(N-k), IGUAL que `newey` de Stata
    # (statsmodels con HAC por defecto usa Normal; Stata usa t).  use_correction=True
    # replica el ajuste de muestra chica del kernel de Bartlett de Stata.
    res = model.fit(
        cov_type="HAC",
        cov_kwds={"maxlags": lag, "use_correction": True},
        use_t=True,
    )
    return res


def ols_plain(y, X, add_const=True):
    """MCO clasico (para diagnosticos que corren tras regress, no tras newey)."""
    yv = pd.Series(y).astype(float)
    Xv = pd.DataFrame(X).astype(float)
    d = pd.concat([yv, Xv], axis=1).dropna()
    yy = d.iloc[:, 0]
    XX = d.iloc[:, 1:]
    if add_const:
        XX = sm.add_constant(XX, has_constant="add")
    return sm.OLS(yy, XX).fit()


def coef_table(res, order=None):
    """Extrae coeficientes, SE, t, p, IC95 de un resultado de statsmodels."""
    ci = res.conf_int(0.05)
    rows = []
    names = order if order is not None else res.params.index
    for nm in names:
        rows.append(
            {
                "variable": nm,
                "coef": res.params[nm],
                "SE": res.bse[nm],
                "t": res.tvalues[nm],
                "p": res.pvalues[nm],
                "ci_low": ci.loc[nm, 0],
                "ci_high": ci.loc[nm, 1],
            }
        )
    return pd.DataFrame(rows)


def one_sided_mean_test(series, lag=HAC_LAG):
    """
    Test t unilateral de la media con SE HAC.
    H0: mu <= 0  vs  HA: mu > 0.  Regresion sobre constante con newey.
    Devuelve media, SE, t, p_unilateral, IC95 bilateral, %>0, n.
    """
    s = pd.Series(series).dropna().astype(float)
    n = len(s)
    X = np.ones(n)
    res = sm.OLS(s.values, X).fit(
        cov_type="HAC", cov_kwds={"maxlags": lag, "use_correction": True},
        use_t=True,
    )
    mean = res.params[0]
    se = res.bse[0]
    tstat = res.tvalues[0]
    # p unilateral para HA: mu>0
    dfree = int(res.df_resid)
    p_one = stats.t.sf(tstat, dfree)  # P(T > t)
    ci = res.conf_int(0.05)
    pct_pos = float((s > 0).mean() * 100.0)
    return {
        "n": n,
        "media": mean,
        "SE_HAC": se,
        "t": tstat,
        "p_unilateral_mu>0": p_one,
        "ci95_low": ci[0][0],
        "ci95_high": ci[0][1],
        "pct_trim_pos": pct_pos,
    }


# ==============================================================================
# ETAPA 1 - Series de retorno + descriptivas + grafico
# ==============================================================================
def etapa1(df):
    print("\n" + "=" * 78)
    print(" ETAPA 1 - Series de retorno y estadistica descriptiva")
    print("=" * 78)

    subsets = {
        "Total (2004-2024)": df,
        "Pre-2011 (2004-2010)": df[df["D2011"] == 0],
        "Post-2011 (2011-2024)": df[df["D2011"] == 1],
    }
    rows = []
    for name, sub in subsets.items():
        for var in ["DR", "DRcons"]:
            d = describe(sub[var], var)
            d["subperiodo"] = name
            rows.append(d)
    tab = pd.DataFrame(rows)[
        ["subperiodo", "serie", "n", "media", "mediana", "sd", "CV",
         "min", "max", "asimetria", "curtosis"]
    ]
    tab.to_csv(os.path.join(TAB, "etapa1_descriptivas.csv"), index=False)
    print(tab.round(4).to_string(index=False))

    # numeros clave para cross-check
    for name, sub in subsets.items():
        key = name.split()[0].lower()
        XCHECK[f"e1_DR_media_{key}"] = float(sub["DR"].dropna().mean())
        XCHECK[f"e1_DR_sd_{key}"] = float(sub["DR"].dropna().std(ddof=1))
        XCHECK[f"e1_DRcons_media_{key}"] = float(sub["DRcons"].dropna().mean())

    # --- Grafico de linea de DR con marca del quiebre 2011 ---
    plot_df = df.dropna(subset=["DR"]).copy()
    x = plot_df["t"].values
    fig, ax = plt.subplots(figsize=(11, 5.2))
    ax.plot(x, plot_df["DR"].values, color="#1f4e79", lw=1.6, marker="o",
            ms=3.5, label="DR (diferencial de retorno)")
    ax.axhline(0, color="#888888", lw=0.9, ls="--")
    tb = q_index(YEAR_BREAK, 1)
    ax.axvline(tb, color="#c00000", lw=1.6, ls="-",
               label="Quiebre regimen 2011 (Metas de Inflacion)")
    # medias por subperiodo
    m_pre = plot_df.loc[plot_df["D2011"] == 0, "DR"].mean()
    m_post = plot_df.loc[plot_df["D2011"] == 1, "DR"].mean()
    ax.hlines(m_pre, x.min(), tb, color="#2e75b6", lw=2.2, alpha=0.9,
              label=f"Media pre-2011 = {m_pre:.2f}")
    ax.hlines(m_post, tb, x.max(), color="#548235", lw=2.2, alpha=0.9,
              label=f"Media post-2011 = {m_post:.2f}")
    # etiquetas de eje x como trimestres
    yrs = list(range(2004, 2025, 2))
    ax.set_xticks([q_index(y, 1) for y in yrs])
    ax.set_xticklabels([str(y) for y in yrs])
    ax.set_xlabel("Trimestre")
    ax.set_ylabel("DR (puntos porcentuales, trimestral)")
    ax.set_title("Etapa 1 - Diferencial de retorno DR (Gs vs USD) con marca 2011")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "etapa1_DR_serie.png"), dpi=130)
    plt.close(fig)

    # grafico comparativo DR vs DRcons
    fig, ax = plt.subplots(figsize=(11, 5.2))
    ax.plot(x, plot_df["DR"].values, color="#1f4e79", lw=1.5, label="DR (estandar)")
    ax.plot(plot_df["t"].values, plot_df["DRcons"].values, color="#c55a11",
            lw=1.5, ls="--", label="DRcons (conservadora)")
    ax.axhline(0, color="#888888", lw=0.9, ls=":")
    ax.axvline(tb, color="#c00000", lw=1.4, label="Quiebre 2011")
    ax.set_xticks([q_index(y, 1) for y in yrs])
    ax.set_xticklabels([str(y) for y in yrs])
    ax.set_xlabel("Trimestre")
    ax.set_ylabel("Puntos porcentuales (trimestral)")
    ax.set_title("Etapa 1 - DR vs DRcons")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "etapa1_DR_vs_DRcons.png"), dpi=130)
    plt.close(fig)

    return tab


# ==============================================================================
# ETAPA 2 - Coeficiente de Fama + validacion conservadora
# ==============================================================================
def etapa2(df):
    print("\n" + "=" * 78)
    print(" ETAPA 2 - Coeficiente de Fama + prueba conservadora (DRcons)")
    print("=" * 78)

    # --- Regresion de Fama:  depLead = alpha + beta*difq + eps  (HAC) ---
    # IMPORTANTE (correccion QA 2026-07-04): el regresor es difq = difn/4, el
    # diferencial TRIMESTRALIZADO, en la misma frecuencia que depLead (dep.
    # trimestral). Solo asi H0: beta=1 es la nula UIP correcta. Regresar sobre
    # difn (anualizado) implica beta_UIP=0.25 y sobre-rechaza beta=1.
    # UIP => beta=1.  beta<1 (o negativo) => falla la paridad descubierta.
    samples = {
        "Total": df,
        "Pre-2011": df[df["D2011"] == 0],
        "Post-2011": df[df["D2011"] == 1],
    }
    rows = []
    for name, sub in samples.items():
        res = newey_ols(sub["depLead"], sub[["difq"]], lag=HAC_LAG)
        beta = res.params["difq"]
        se = res.bse["difq"]
        # test H0: beta = 1
        from statsmodels.stats.contrast import ContrastResults  # noqa
        t_b1 = (beta - 1.0) / se
        p_b1 = 2 * stats.t.sf(abs(t_b1), int(res.df_resid))
        rows.append(
            {
                "muestra": name,
                "n": int(res.nobs),
                "alpha": res.params["const"],
                "beta": beta,
                "SE_beta": se,
                "p_beta=0": res.pvalues["difq"],
                "t_(beta=1)": t_b1,
                "p_(beta=1)": p_b1,
                "R2": res.rsquared,
            }
        )
        XCHECK[f"e2_fama_beta_{name.lower()}"] = float(beta)
    fama = pd.DataFrame(rows)
    fama.to_csv(os.path.join(TAB, "etapa2_fama.csv"), index=False)
    print("\n-- Regresion de Fama: depLead = alpha + beta*difq (HAC Newey-West) --")
    print("   UIP <=> beta=1.  beta<1 => falla la paridad descubierta.")
    print(fama.round(4).to_string(index=False))

    # --- Prueba conservadora: H0 mu(DRcons)<=0 vs HA>0 (HAC) ---
    print("\n-- Prueba unilateral H0: mu(DRcons)<=0  vs  HA: mu>0  (HAC) --")
    rows2 = []
    for name, sub in samples.items():
        r = one_sided_mean_test(sub["DRcons"], lag=HAC_LAG)
        r["muestra"] = name
        rows2.append(r)
        XCHECK[f"e2_DRcons_media_{name.lower()}"] = float(r["media"])
        XCHECK[f"e2_DRcons_pct_pos_{name.lower()}"] = float(r["pct_trim_pos"])
    cons = pd.DataFrame(rows2)[
        ["muestra", "n", "media", "SE_HAC", "t", "p_unilateral_mu>0",
         "ci95_low", "ci95_high", "pct_trim_pos"]
    ]
    cons.to_csv(os.path.join(TAB, "etapa2_drcons_test.csv"), index=False)
    print(cons.round(4).to_string(index=False))

    # --- Comparacion DR vs DRcons (medias) ---
    cmp_rows = []
    for name, sub in samples.items():
        cmp_rows.append(
            {
                "muestra": name,
                "media_DR": sub["DR"].dropna().mean(),
                "media_DRcons": sub["DRcons"].dropna().mean(),
                "pct_DR_pos": float((sub["DR"].dropna() > 0).mean() * 100),
                "pct_DRcons_pos": float((sub["DRcons"].dropna() > 0).mean() * 100),
            }
        )
    cmp = pd.DataFrame(cmp_rows)
    cmp.to_csv(os.path.join(TAB, "etapa2_DR_vs_DRcons.csv"), index=False)
    print("\n-- Comparacion DR vs DRcons --")
    print(cmp.round(4).to_string(index=False))

    # --- Robustez: Fama con series fin-de-trimestre (difq_fdt) ---
    rows3 = []
    for name, sub in samples.items():
        res = newey_ols(sub["depLead"], sub[["difq_fdt"]], lag=HAC_LAG)
        rows3.append(
            {"muestra": name, "n": int(res.nobs),
             "beta_fdt": res.params["difq_fdt"], "SE": res.bse["difq_fdt"],
             "p": res.pvalues["difq_fdt"]}
        )
    fama_fdt = pd.DataFrame(rows3)
    fama_fdt.to_csv(os.path.join(TAB, "etapa2_fama_fdt_robustez.csv"), index=False)
    print("\n-- Robustez Fama (series fin-de-trimestre i_*_fdt) --")
    print(fama_fdt.round(4).to_string(index=False))

    return fama, cons, cmp


# ==============================================================================
# ETAPA 3 - Estacionariedad + quiebre estructural
# ==============================================================================
def bic_lag_adf(series, maxlag=8):
    """Elige lags del ADF por BIC (equivalente al criterio autolag='BIC')."""
    s = series.dropna().astype(float)
    best = adfuller(s, maxlag=maxlag, regression="c", autolag="BIC")
    return best


def phillips_perron(s, trend="c"):
    """
    Test de Phillips-Perron (estadistico Z_tau) con constante.
    Replica `pperron` de Stata:
      1) DF sin aumentar:  Dy_t = a + rho*y_{t-1} + u_t
      2) t_rho crudo + correccion no parametrica con varianza de largo plazo
         estimada por kernel de Bartlett (lags Newey-West = floor(4*(T/100)^{2/9})).
    """
    s = np.asarray(s, float)
    y = s[1:]
    ylag = s[:-1]
    dy = np.diff(s)
    n = len(dy)
    if trend == "c":
        X = np.column_stack([np.ones(n), ylag])
    else:  # 'ct'
        X = np.column_stack([np.ones(n), np.arange(1, n + 1), ylag])
    beta = np.linalg.lstsq(X, dy, rcond=None)[0]
    resid = dy - X @ beta
    XtX_inv = np.linalg.inv(X.T @ X)
    rho_idx = 1 if trend == "c" else 2
    s2 = (resid @ resid) / (n - X.shape[1])
    se_rho = np.sqrt(s2 * XtX_inv[rho_idx, rho_idx])
    t_rho = beta[rho_idx] / se_rho

    # varianza de corto plazo (sigma^2) y de largo plazo (lambda^2) via Bartlett
    l = int(np.floor(4 * (n / 100.0) ** (2.0 / 9.0)))
    gamma0 = (resid @ resid) / n
    lam2 = gamma0
    for j in range(1, l + 1):
        w = 1.0 - j / (l + 1.0)
        gj = (resid[j:] @ resid[:-j]) / n
        lam2 += 2.0 * w * gj

    # correccion de Phillips-Perron sobre el estadistico t
    #   Z_t = sqrt(gamma0/lam2)*t_rho - (lam2-gamma0)/(2*lam2) * (n*se_rho/sqrt(s2))
    term2 = (lam2 - gamma0) / (2.0 * lam2) * (n * se_rho / np.sqrt(gamma0))
    z_tau = np.sqrt(gamma0 / lam2) * t_rho - term2

    # p-valor via superficie de MacKinnon (usar tabla de statsmodels ADF)
    from statsmodels.tsa.adfvalues import mackinnonp

    p = mackinnonp(z_tau, regression=trend, N=1)
    return float(z_tau), float(p), l


def etapa3(df):
    print("\n" + "=" * 78)
    print(" ETAPA 3 - Estacionariedad (ADF/PP) + quiebre estructural")
    print("=" * 78)

    s = df["DR"].dropna().astype(float).reset_index(drop=True)

    # --- ADF con constante (lags por BIC) ---
    adf_c = adfuller(s, maxlag=8, regression="c", autolag="BIC")
    # --- ADF con constante + tendencia ---
    adf_ct = adfuller(s, maxlag=8, regression="ct", autolag="BIC")
    # --- Phillips-Perron (Z_tau, con constante) ---
    # statsmodels no trae PP; se implementa la formula estandar:
    #   regresion DF sin aumentar:  Dy_t = a + rho*y_{t-1} + u_t
    #   correccion no parametrica de la varianza de largo plazo (Newey-West Bartlett).
    # Reproduce `pperron DR` de Stata (lags Newey-West por defecto ~ 4*(T/100)^{2/9}).
    pp_stat, pp_p, pp_lags = phillips_perron(s)

    adf_rows = [
        {"test": "ADF (const)", "stat": adf_c[0], "p": adf_c[1],
         "lags": adf_c[2], "n": adf_c[3],
         "cv_5pct": adf_c[4]["5%"]},
        {"test": "ADF (const+trend)", "stat": adf_ct[0], "p": adf_ct[1],
         "lags": adf_ct[2], "n": adf_ct[3],
         "cv_5pct": adf_ct[4]["5%"]},
        {"test": "Phillips-Perron (const)", "stat": pp_stat,
         "p": pp_p, "lags": pp_lags, "n": len(s) - 1,
         "cv_5pct": adf_c[4]["5%"]},
    ]
    adf_tab = pd.DataFrame(adf_rows)
    adf_tab.to_csv(os.path.join(TAB, "etapa3_raiz_unitaria.csv"), index=False)
    print("\n-- Raiz unitaria sobre DR (H0: raiz unitaria / no estacionaria) --")
    print(adf_tab.round(4).to_string(index=False))
    XCHECK["e3_adf_c_stat"] = float(adf_c[0])
    XCHECK["e3_adf_c_p"] = float(adf_c[1])
    XCHECK["e3_adf_c_lags"] = int(adf_c[2])
    XCHECK["e3_pp_stat"] = float(pp_stat)

    # --- Test t de diferencia de medias pre/post 2011 (HAC) ---
    # regresion DR = a + b*D2011 ; b = diferencia de medias post - pre
    dd = df.dropna(subset=["DR"])
    res = newey_ols(dd["DR"], dd[["D2011"]], lag=HAC_LAG)
    diff_rows = [
        {"parametro": "media pre-2011 (const)", "valor": res.params["const"],
         "SE": res.bse["const"], "t": res.tvalues["const"], "p": res.pvalues["const"]},
        {"parametro": "diferencia post-pre (D2011)", "valor": res.params["D2011"],
         "SE": res.bse["D2011"], "t": res.tvalues["D2011"], "p": res.pvalues["D2011"]},
    ]
    diff_tab = pd.DataFrame(diff_rows)
    diff_tab.to_csv(os.path.join(TAB, "etapa3_dif_medias.csv"), index=False)
    print("\n-- Diferencia de medias DR pre/post 2011 (HAC) --")
    print(diff_tab.round(4).to_string(index=False))
    XCHECK["e3_dif_medias_D2011"] = float(res.params["D2011"])
    XCHECK["e3_dif_medias_p"] = float(res.pvalues["D2011"])

    # --- Quiebre estructural endogeno: Quandt-Andrews (sup-Wald) ---
    # Se replica lo que hace `estat sbsingle`: barre fechas de quiebre en el
    # 15%-85% central y toma el max estadistico de Wald (sup-Wald / sup-F).
    # Modelo base: DR = c + phi*L.DR  (mismo que el do-file: regress DR L.DR).
    dr_lag = pd.DataFrame({"DR": df["DR"], "L1": df["DR"].shift(1),
                           "t": df["t"], "anio": df["anio"],
                           "trim": df["trimestre"]}).dropna().reset_index(drop=True)
    yv = dr_lag["DR"].values
    Xbase = sm.add_constant(dr_lag[["L1"]].values)
    n = len(yv)
    trim_frac = 0.15
    lo = int(np.floor(n * trim_frac))
    hi = int(np.ceil(n * (1 - trim_frac)))
    sup_wald = -np.inf
    sup_idx = None
    k = Xbase.shape[1]
    wald_series = []
    for i in range(lo, hi):
        d = np.zeros(n)
        d[i:] = 1.0
        Xi = np.column_stack([Xbase, d, d * dr_lag["L1"].values])
        try:
            r = sm.OLS(yv, Xi).fit()
            # Wald de que los 2 coef de interaccion del quiebre son 0
            Rmat = np.zeros((2, Xi.shape[1]))
            Rmat[0, k] = 1.0
            Rmat[1, k + 1] = 1.0
            wt = r.wald_test(Rmat, use_f=False, scalar=True)
            wstat = float(np.asarray(wt.statistic).ravel()[0])
        except Exception:
            wstat = np.nan
        wald_series.append((dr_lag["anio"][i], dr_lag["trim"][i], wstat))
        if np.isfinite(wstat) and wstat > sup_wald:
            sup_wald = wstat
            sup_idx = i
    brk_year = int(dr_lag["anio"][sup_idx])
    brk_trim = int(dr_lag["trim"][sup_idx])
    # p-valor aproximado de Andrews(1993)/Hansen(1997) via distribucion asintotica.
    # Se reporta el estadistico y la fecha; el p exacto lo entrega Stata sbsingle.
    print("\n-- Quiebre unico endogeno (Quandt-Andrews sup-Wald sobre DR=c+phi*L.DR) --")
    print(f"   sup-Wald = {sup_wald:.4f}  en  {brk_year}Q{brk_trim}")
    XCHECK["e3_supwald_stat"] = float(sup_wald)
    XCHECK["e3_supwald_fecha"] = f"{brk_year}Q{brk_trim}"

    # tabla de la trayectoria de Wald (para inspeccion)
    wtab = pd.DataFrame(wald_series, columns=["anio", "trimestre", "wald"])
    wtab.to_csv(os.path.join(TAB, "etapa3_supwald_trayectoria.csv"), index=False)

    # --- Bai-Perron (quiebres multiples) via ruptures sobre el nivel de DR ---
    bp_rows = []
    if HAS_RUPTURES:
        sig = s.values.reshape(-1, 1)
        # Dynp con modelo de media (l2); minimo 5 obs por segmento (~15% de 81)
        fechas = df.dropna(subset=["DR"])[["anio", "trimestre"]].reset_index(drop=True)
        for nbk in (1, 2, 3):
            try:
                algo = rpt.Dynp(model="l2", min_size=5, jump=1).fit(sig)
                bkps = algo.predict(n_bkps=nbk)
                # bkps incluye el largo total al final; los quiebres son los internos
                internos = [b for b in bkps if b < len(sig)]
                puntos = []
                for b in internos:
                    yb = int(fechas["anio"][b])
                    tb_ = int(fechas["trimestre"][b])
                    puntos.append(f"{yb}Q{tb_}")
                bp_rows.append({"n_quiebres": nbk, "fechas": "; ".join(puntos)})
            except Exception as e:
                bp_rows.append({"n_quiebres": nbk, "fechas": f"error: {e}"})
        # seleccion por penalizacion (Pelt/BIC-like) para numero endogeno de quiebres
        try:
            algo_pelt = rpt.Pelt(model="l2", min_size=5, jump=1).fit(sig)
            # penalizacion ~ log(n)*sigma^2 (BIC)
            sigma2 = np.var(sig)
            pen = np.log(len(sig)) * sigma2
            bkps_pelt = algo_pelt.predict(pen=pen)
            internos = [b for b in bkps_pelt if b < len(sig)]
            fechas_pelt = []
            for b in internos:
                fechas_pelt.append(f"{int(fechas['anio'][b])}Q{int(fechas['trimestre'][b])}")
            bp_rows.append({"n_quiebres": f"PELT/BIC (endogeno={len(internos)})",
                            "fechas": "; ".join(fechas_pelt) if fechas_pelt else "ninguno"})
            XCHECK["e3_bp_pelt_nquiebres"] = len(internos)
            XCHECK["e3_bp_pelt_fechas"] = "; ".join(fechas_pelt)
        except Exception:
            pass
    bp_tab = pd.DataFrame(bp_rows)
    bp_tab.to_csv(os.path.join(TAB, "etapa3_baiperron_ruptures.csv"), index=False)
    print("\n-- Bai-Perron (ruptures) sobre nivel de DR --")
    print(bp_tab.to_string(index=False))

    return adf_tab, diff_tab, bp_tab


# ==============================================================================
# ETAPA 4 - Descomposicion de Fisher
# ==============================================================================
def etapa4(df):
    print("\n" + "=" * 78)
    print(" ETAPA 4 - Descomposicion de Fisher (prima inflacion vs riesgo)")
    print("=" * 78)

    # difn = InflDiff + rp   (por construccion, rp = difn - InflDiff)
    d = df.dropna(subset=["difn", "InflDiff", "rp"]).copy()

    rows = []
    for name, sub in {"Pre-2011": d[d["D2011"] == 0],
                      "Post-2011": d[d["D2011"] == 1],
                      "Total": d}.items():
        difn_m = sub["difn"].mean()
        infl_m = sub["InflDiff"].mean()
        rp_m = sub["rp"].mean()
        rows.append(
            {
                "subperiodo": name,
                "n": len(sub),
                "difn_medio": difn_m,
                "InflDiff_medio": infl_m,
                "rp_medio": rp_m,
                "prop_inflacion_%": 100 * infl_m / difn_m if difn_m != 0 else np.nan,
                "prop_riesgo_%": 100 * rp_m / difn_m if difn_m != 0 else np.nan,
                "sd_InflDiff": sub["InflDiff"].std(ddof=1),
                "sd_rp": sub["rp"].std(ddof=1),
            }
        )
    fisher = pd.DataFrame(rows)
    fisher.to_csv(os.path.join(TAB, "etapa4_fisher.csv"), index=False)
    print("\n-- Descomposicion difn = InflDiff (prima inflacion) + rp (prima riesgo) --")
    print(fisher.round(4).to_string(index=False))

    for _, r in fisher.iterrows():
        k = r["subperiodo"].split("-")[0].lower() if "-" in r["subperiodo"] else r["subperiodo"].lower()
        XCHECK[f"e4_prop_infl_{k}"] = float(r["prop_inflacion_%"])
        XCHECK[f"e4_prop_riesgo_{k}"] = float(r["prop_riesgo_%"])

    # --- prueba t de diferencia de medias pre/post (HAC) sobre rp e InflDiff ---
    rows2 = []
    for var in ["rp", "InflDiff", "difn"]:
        res = newey_ols(d[var], d[["D2011"]], lag=HAC_LAG)
        rows2.append(
            {
                "variable": var,
                "media_pre (const)": res.params["const"],
                "dif_post-pre (D2011)": res.params["D2011"],
                "SE": res.bse["D2011"],
                "t": res.tvalues["D2011"],
                "p": res.pvalues["D2011"],
            }
        )
    fisher_t = pd.DataFrame(rows2)
    fisher_t.to_csv(os.path.join(TAB, "etapa4_fisher_ttest.csv"), index=False)
    print("\n-- Prueba t de diferencia de medias pre/post 2011 (HAC) --")
    print(fisher_t.round(4).to_string(index=False))
    XCHECK["e4_rp_dif_D2011"] = float(fisher_t.loc[fisher_t["variable"] == "rp", "dif_post-pre (D2011)"].values[0])
    XCHECK["e4_rp_dif_p"] = float(fisher_t.loc[fisher_t["variable"] == "rp", "p"].values[0])

    # --- grafico apilado de la descomposicion por subperiodo ---
    fig, ax = plt.subplots(figsize=(7.5, 5))
    cats = ["Pre-2011", "Post-2011"]
    infl_vals = [fisher.loc[fisher["subperiodo"] == c, "InflDiff_medio"].values[0] for c in cats]
    rp_vals = [fisher.loc[fisher["subperiodo"] == c, "rp_medio"].values[0] for c in cats]
    x = np.arange(len(cats))
    ax.bar(x, infl_vals, 0.55, label="Prima de inflacion (InflDiff)", color="#2e75b6")
    ax.bar(x, rp_vals, 0.55, bottom=infl_vals, label="Prima de riesgo (rp)", color="#c55a11")
    for i, c in enumerate(cats):
        tot = infl_vals[i] + rp_vals[i]
        ax.text(i, tot + 0.1, f"difn={tot:.2f}", ha="center", fontsize=9, fontweight="bold")
        ax.text(i, infl_vals[i] / 2, f"{100*infl_vals[i]/tot:.0f}%", ha="center",
                color="white", fontsize=9)
        ax.text(i, infl_vals[i] + rp_vals[i] / 2, f"{100*rp_vals[i]/tot:.0f}%", ha="center",
                color="white", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(cats)
    ax.set_ylabel("Puntos porcentuales (% anual)")
    ax.set_title("Etapa 4 - Descomposicion de Fisher del diferencial nominal")
    ax.legend(fontsize=9)
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "etapa4_fisher_descomposicion.png"), dpi=130)
    plt.close(fig)

    return fisher, fisher_t


# ==============================================================================
# ETAPA 5 - Modelo MCO-HAC con interacciones (NUCLEO del contraste)
# ==============================================================================
REGRESSORS = ["InflDiff", "VolTC", "Iliq", "D2011", "D_VolTC", "D_Iliq"]
COEF_LABELS = {
    "const": "beta0 (const)",
    "InflDiff": "beta1 InflDiff",
    "VolTC": "beta2 VolTC",
    "Iliq": "beta3 Iliq",
    "D2011": "beta4 D2011",
    "D_VolTC": "beta5 D2011xVolTC",
    "D_Iliq": "beta6 D2011xIliq",
}


def etapa5(df):
    print("\n" + "=" * 78)
    print(" ETAPA 5 - MCO-HAC con interacciones (NUCLEO)")
    print("=" * 78)
    print(" DR = b0 + b1 InflDiff + b2 VolTC + b3 Iliq + b4 D2011")
    print("       + b5 (D2011xVolTC) + b6 (D2011xIliq) + e     [HAC lag=%d]" % HAC_LAG)
    print(" Hipotesis de recomposicion:  HA b5<0 (VolTC pierde peso) y b6>0 (Iliq gana peso)")

    results = {}
    for dep in ["DR", "DRcons"]:
        res = newey_ols(df[dep], df[REGRESSORS], lag=HAC_LAG)
        order = ["const"] + REGRESSORS
        ct = coef_table(res, order=order)
        ct["variable"] = ct["variable"].map(COEF_LABELS)
        # p unilateral para b5 (HA<0) y b6 (HA>0)
        dfree = int(res.df_resid)
        p1_b5 = stats.t.cdf(res.tvalues["D_VolTC"], dfree)   # HA: b5<0 -> P(T<t)
        p1_b6 = stats.t.sf(res.tvalues["D_Iliq"], dfree)     # HA: b6>0 -> P(T>t)
        ct["p_unilateral_HA"] = np.nan
        ct.loc[ct["variable"] == COEF_LABELS["D_VolTC"], "p_unilateral_HA"] = p1_b5
        ct.loc[ct["variable"] == COEF_LABELS["D_Iliq"], "p_unilateral_HA"] = p1_b6
        ct.to_csv(os.path.join(TAB, f"etapa5_modelo_{dep}.csv"), index=False)
        print(f"\n-- Modelo con dependiente = {dep}  (n={int(res.nobs)}, R2={res.rsquared:.3f}) --")
        print(ct.round(4).to_string(index=False))
        print(f"   Contraste b5 (D2011xVolTC): coef={res.params['D_VolTC']:.4f}, "
              f"p_2c={res.pvalues['D_VolTC']:.4f}, p_1c(HA<0)={p1_b5:.4f}  "
              f"[esperado <0]")
        print(f"   Contraste b6 (D2011xIliq):  coef={res.params['D_Iliq']:.4f}, "
              f"p_2c={res.pvalues['D_Iliq']:.4f}, p_1c(HA>0)={p1_b6:.4f}  "
              f"[esperado >0]")
        results[dep] = res
        # cross-check
        for v in REGRESSORS:
            XCHECK[f"e5_{dep}_b_{v}"] = float(res.params[v])
            XCHECK[f"e5_{dep}_p_{v}"] = float(res.pvalues[v])
        XCHECK[f"e5_{dep}_n"] = int(res.nobs)
        XCHECK[f"e5_{dep}_r2"] = float(res.rsquared)

    return results


# ==============================================================================
# ETAPA 6 - Diagnostico + robustez
# ==============================================================================
def etapa6(df):
    print("\n" + "=" * 78)
    print(" ETAPA 6 - Diagnostico + robustez")
    print("=" * 78)

    # --- OLS clasico paralelo (para diagnosticos que corren tras regress) ---
    res_ols = ols_plain(df["DR"], df[REGRESSORS])

    # Breusch-Godfrey (autocorrelacion) hasta 4 rezagos
    diag_rows = []
    for L in (1, 2, 3, 4):
        bg = acorr_breusch_godfrey(res_ols, nlags=L)
        diag_rows.append({"test": f"Breusch-Godfrey (lag {L})",
                          "stat_LM": bg[0], "p": bg[1]})
    # White (heterocedasticidad)
    wh = het_white(res_ols.resid, res_ols.model.exog)
    diag_rows.append({"test": "White (heterocedasticidad)", "stat_LM": wh[0], "p": wh[1]})
    # Ramsey RESET (forma funcional)
    try:
        reset = linear_reset(res_ols, power=2, use_f=True)
        diag_rows.append({"test": "Ramsey RESET (power 2, F)",
                          "stat_LM": float(reset.statistic), "p": float(reset.pvalue)})
        reset3 = linear_reset(res_ols, power=3, use_f=True)
        diag_rows.append({"test": "Ramsey RESET (power 2-3, F)",
                          "stat_LM": float(reset3.statistic), "p": float(reset3.pvalue)})
    except Exception as e:
        diag_rows.append({"test": "Ramsey RESET", "stat_LM": np.nan, "p": np.nan})
    diag = pd.DataFrame(diag_rows)
    diag.to_csv(os.path.join(TAB, "etapa6_diagnosticos.csv"), index=False)
    print("\n-- Diagnosticos sobre el OLS clasico (DR) --")
    print(diag.round(4).to_string(index=False))
    XCHECK["e6_bg4_p"] = float(diag.loc[diag["test"] == "Breusch-Godfrey (lag 4)", "p"].values[0])
    XCHECK["e6_white_p"] = float(diag.loc[diag["test"] == "White (heterocedasticidad)", "p"].values[0])
    XCHECK["e6_reset_p"] = float(diag.loc[diag["test"] == "Ramsey RESET (power 2, F)", "p"].values[0])

    # --- VIF / colinealidad entre regresores (limitacion declarada) ---
    Xv = df[REGRESSORS].dropna()
    Xc = sm.add_constant(Xv)
    vif_rows = []
    for i, nm in enumerate(Xc.columns):
        if nm == "const":
            continue
        vif_rows.append({"variable": nm, "VIF": variance_inflation_factor(Xc.values, i)})
    vif = pd.DataFrame(vif_rows)
    vif.to_csv(os.path.join(TAB, "etapa6_vif.csv"), index=False)
    print("\n-- VIF (multicolinealidad; alta = limitacion real del diseno T=84) --")
    print(vif.round(3).to_string(index=False))
    # correlaciones D2011 ~ D_Iliq ~ D_VolTC
    corr = df[["D2011", "D_Iliq", "D_VolTC", "Iliq"]].corr()
    corr.to_csv(os.path.join(TAB, "etapa6_correlaciones.csv"))
    print("\n-- Correlaciones D2011 / D_Iliq / D_VolTC / Iliq --")
    print(corr.round(3).to_string())
    XCHECK["e6_vif_max"] = float(vif["VIF"].max())
    XCHECK["e6_corr_D2011_DIliq"] = float(corr.loc["D2011", "D_Iliq"])

    # --- Robustez (ii): DRcons como dependiente (ya en Etapa 5, se re-lista aqui) ---
    res_cons = newey_ols(df["DRcons"], df[REGRESSORS], lag=HAC_LAG)

    # --- Robustez (iii): submuestra 2011-2024 (post homogeneo) ---
    #     En post no hay variacion de D2011 -> se estima el modelo SIN dummies:
    #     DR = a + b1 InflDiff + b2 VolTC + b3 Iliq
    post = df[df["D2011"] == 1]
    base_reg = ["InflDiff", "VolTC", "Iliq"]
    res_post = newey_ols(post["DR"], post[base_reg], lag=HAC_LAG)
    res_post_cons = newey_ols(post["DRcons"], post[base_reg], lag=HAC_LAG)

    rob_rows = []
    for lbl, res, regs in [
        ("Full DR (6 regs)", None, None),
        ("Full DRcons (6 regs)", res_cons, REGRESSORS),
        ("Post-2011 DR (3 regs)", res_post, base_reg),
        ("Post-2011 DRcons (3 regs)", res_post_cons, base_reg),
    ]:
        if res is None:
            continue
        d = {"modelo": lbl, "n": int(res.nobs), "R2": round(res.rsquared, 3)}
        for v in ["InflDiff", "VolTC", "Iliq", "D_VolTC", "D_Iliq"]:
            if v in res.params.index:
                d[f"b_{v}"] = round(res.params[v], 4)
                d[f"p_{v}"] = round(res.pvalues[v], 4)
            else:
                d[f"b_{v}"] = np.nan
                d[f"p_{v}"] = np.nan
        rob_rows.append(d)
    rob = pd.DataFrame(rob_rows)
    rob.to_csv(os.path.join(TAB, "etapa6_robustez.csv"), index=False)
    print("\n-- Robustez: DRcons y submuestra post-2011 (sin dummies) --")
    print(rob.to_string(index=False))

    # --- Robustez (i): brecha de encaje -> NO FACTIBLE ---
    nota_encaje = (
        "Robustez (i) brecha de encaje Gs/USD: NO FACTIBLE. Las variables "
        "encaje_Gs / encaje_USD NO estan en la base porque el recopilador no las "
        "obtuvo de fuente oficial (BCP). Se declara como limitacion. dolariz (Iliq) "
        "se mantiene como proxy primaria de iliquidez, que es lo que la tesis define "
        "como principal."
    )
    with open(os.path.join(TAB, "etapa6_encaje_NO_FACTIBLE.txt"), "w", encoding="utf-8") as f:
        f.write(nota_encaje + "\n")
    print("\n-- Robustez (i) brecha de encaje --")
    print("   " + nota_encaje)

    return diag, vif, rob


# ==============================================================================
# MAIN
# ==============================================================================
def main():
    print("#" * 78)
    print("# ANALISIS ECONOMETRICO - TESIS CDA Gs vs USD (Python == Stata)")
    print("# HAC Newey-West lag =", HAC_LAG, " | quiebre regimen =", YEAR_BREAK)
    print("#" * 78)

    df = load_and_build()
    print(f"\nBase cargada: {df.shape[0]} filas (2004Q1-2024Q4).")
    print(f"Missing i_Gs/i_USD (exclusion listwise, sin imputar): "
          f"{int(df['i_Gs'].isna().sum())} -> "
          f"{df.loc[df['i_Gs'].isna(), 'fecha'].tolist()}")

    etapa1(df)
    etapa2(df)
    etapa3(df)
    etapa4(df)
    etapa5(df)
    etapa6(df)

    # volcar numeros clave para el cross-check con Stata
    with open(os.path.join(OUT, "python_xcheck.json"), "w", encoding="utf-8") as f:
        json.dump(XCHECK, f, indent=2, ensure_ascii=False)
    print("\n" + "=" * 78)
    print(" Numeros clave para cross-check -> output/python_xcheck.json")
    print(" Tablas -> output/tablas/*.csv   |   Figuras -> output/figuras/*.png")
    print("=" * 78)

    # tambien guardar la base derivada para que Stata pueda re-usar exactamente
    keep = ["anio", "trimestre", "fecha", "t", "i_Gs", "i_USD", "TC", "dolariz",
            "difn", "difq", "dep", "depLead", "DR", "DRcons", "D2011",
            "infl_py", "infl_us", "InflDiff", "VolTC", "rp", "Iliq",
            "D_VolTC", "D_Iliq"]
    df[keep].to_csv(os.path.join(OUT, "base_derivada_python.csv"), index=False)
    print(" Base derivada -> output/base_derivada_python.csv")


if __name__ == "__main__":
    main()
