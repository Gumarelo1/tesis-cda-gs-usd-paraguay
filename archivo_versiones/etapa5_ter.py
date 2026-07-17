# -*- coding: utf-8 -*-
"""
ETAPA 5-TER — Modelo principal de la tesis (Tabla 8.3) y verificaciones.
Productor trazable de:
  - output/tablas/etapa5_modelo_rp_centrado.csv  (Tabla 8.3 del documento)
  - output/tablas/etapa5_revisada_tutor.json     (todas las verificaciones)
Equivalente Python de la ETAPA 5-TER de stata/tesis_cda.do.
"""
import json

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.stattools import adfuller, kpss
import warnings
warnings.filterwarnings("ignore")

BASE = "data/processed/base_cda.csv"
HAC = 4


def sc(v):
    return float(np.asarray(v).squeeze())


def hac(y, X):
    return sm.OLS(y, sm.add_constant(X)).fit(
        cov_type="HAC", cov_kwds={"maxlags": HAC}, use_t=True)


def main() -> None:
    df = pd.read_csv(BASE).sort_values(["anio", "trimestre"]).reset_index(drop=True)
    df["difn"] = df["i_Gs"] - df["i_USD"]
    df["dep"] = df["TC"].pct_change() * 100
    df["infl_py"] = (df["IPC_PY"] / df["IPC_PY"].shift(4) - 1) * 100
    df["infl_us"] = (df["IPC_US"] / df["IPC_US"].shift(4) - 1) * 100
    df["InflDiff"] = df["infl_py"] - df["infl_us"]
    df["rp"] = df["difn"] - df["InflDiff"]
    df["VolTC"] = df["dep"].rolling(8, min_periods=2).std()
    df["Iliq"] = df["dolariz"]
    df["D2011"] = (df["anio"] >= 2011).astype(float)
    df["L_Iliq"] = df["Iliq"].shift(1)

    d = df.dropna(subset=["rp", "VolTC", "Iliq"]).reset_index(drop=True)
    mV, sV = d.VolTC.mean(), d.VolTC.std(ddof=1)
    mI, sI = d.Iliq.mean(), d.Iliq.std(ddof=1)
    d["VolTC_c"] = d.VolTC - mV
    d["Iliq_c"] = d.Iliq - mI
    d["DxV"] = d.D2011 * d.VolTC_c
    d["DxI"] = d.D2011 * d.Iliq_c

    # (1) modelo principal (Tabla 8.3)
    m = hac(d.rp, d[["VolTC_c", "Iliq_c", "D2011", "DxV", "DxI"]])
    w = m.wald_test("DxV = 0, DxI = 0", use_f=True, scalar=True)
    wf = m.wald_test("D2011 = 0, DxV = 0, DxI = 0", use_f=True, scalar=True)
    tabla = pd.DataFrame({
        "coef": ["b0", "b2_VolTC_c", "b3_Iliq_c", "b4_D2011", "b5_DxVolTC", "b6_DxIliq"],
        "estimado": [m.params[v] for v in ["const", "VolTC_c", "Iliq_c", "D2011", "DxV", "DxI"]],
        "SE_HAC": [m.bse[v] for v in ["const", "VolTC_c", "Iliq_c", "D2011", "DxV", "DxI"]],
        "p": [m.pvalues[v] for v in ["const", "VolTC_c", "Iliq_c", "D2011", "DxV", "DxI"]],
    })
    tabla.loc[len(tabla)] = ["Wald_F(b5=b6=0)", sc(w.statistic), np.nan, sc(w.pvalue)]
    tabla.loc[len(tabla)] = ["Chow_F", sc(wf.statistic), np.nan, sc(wf.pvalue)]
    tabla.loc[len(tabla)] = ["n / R2", m.nobs, np.nan, m.rsquared]
    tabla.to_csv("output/tablas/etapa5_modelo_rp_centrado.csv", index=False)
    print(f"(1) Tabla 8.3: D2011={m.params['D2011']:+.3f} (p={m.pvalues['D2011']:.4f}) "
          f"b6={m.params['DxI']:+.4f} (p={m.pvalues['DxI']:.4f}) "
          f"Wald F={sc(w.statistic):.3f} p={sc(w.pvalue):.4f}")

    # (2) dominancia estandarizada
    d["zV"] = (d.VolTC - mV) / sV
    d["zI"] = (d.Iliq - mI) / sI
    d["DzV"] = d.D2011 * d.zV
    d["DzI"] = d.D2011 * d.zI
    mz = hac(d.rp, d[["zV", "zI", "D2011", "DzV", "DzI"]])
    incV = mz.params["zV"] + mz.params["DzV"]
    incI = mz.params["zI"] + mz.params["DzI"]
    tdif = mz.t_test("zI + DzI - zV - DzV = 0")
    print(f"(2) dominancia: incV={incV:+.3f}sd incI={incI:+.3f}sd "
          f"t={sc(tdif.tvalue):+.3f} p={sc(tdif.pvalue):.4f}")

    # (3) integracion (misma muestra que el modelo)
    integ = []
    for v in ["rp", "VolTC", "Iliq", "InflDiff"]:
        s = d[v].dropna()
        a = adfuller(s, autolag="BIC")
        k = kpss(s, regression="c", nlags="auto")
        integ.append({"serie": v, "ADF": a[0], "pADF": a[1], "KPSS": k[0], "pKPSS": k[1]})
    print("(3) integracion: dolarizacion I(1) (limite declarado)")

    # (4) JB + ADF de residuos (descarta espuriedad)
    jb, jbp, _, _ = sm.stats.jarque_bera(m.resid)
    ar = adfuller(m.resid, autolag="BIC")
    print(f"(4) JB={jb:.2f} (p={jbp:.4f}) | ADF residuos={ar[0]:.3f} (p={ar[1]:.5f})")

    # (5) robustez con dolarizacion rezagada
    dl = df.dropna(subset=["rp", "VolTC", "L_Iliq"]).copy()
    dl["VolTC_c"] = dl.VolTC - dl.VolTC.mean()
    dl["LIc"] = dl.L_Iliq - dl.L_Iliq.mean()
    dl["DxV"] = dl.D2011 * dl.VolTC_c
    dl["DxLI"] = dl.D2011 * dl.LIc
    ml = hac(dl.rp, dl[["VolTC_c", "LIc", "D2011", "DxV", "DxLI"]])
    wl = ml.wald_test("DxV = 0, DxLI = 0", use_f=True, scalar=True)
    print(f"(5) rezago: Wald F={sc(wl.statistic):.3f} p={sc(wl.pvalue):.4f}")

    # (6) regimen gradual
    def inten(r):
        if r.anio < 2011:
            return 0.0
        if r.anio >= 2017:
            return 1.0
        return ((r.anio - 2011) * 4 + (r.trimestre - 1) + 1) / 24
    d["Reg"] = d.apply(inten, axis=1)
    d["RxV"] = d.Reg * d.VolTC_c
    d["RxI"] = d.Reg * d.Iliq_c
    mg = hac(d.rp, d[["VolTC_c", "Iliq_c", "Reg", "RxV", "RxI"]])
    wg = mg.wald_test("RxV = 0, RxI = 0", use_f=True, scalar=True)
    print(f"(6) gradual: Wald F={sc(wg.statistic):.3f} p={sc(wg.pvalue):.4f}")

    out = {
        "modelo_centrado": {v: [float(m.params[v]), float(m.bse[v]), float(m.pvalues[v])]
                            for v in m.params.index},
        "R2": float(m.rsquared), "n": int(m.nobs),
        "wald_F": sc(w.statistic), "wald_p": sc(w.pvalue),
        "chow_F": sc(wf.statistic), "chow_p": sc(wf.pvalue),
        "dominancia": {"incV_post_sd": float(incV), "incI_post_sd": float(incI),
                       "t_dif": sc(tdif.tvalue), "p_dif": sc(tdif.pvalue)},
        "integracion": integ,
        "jb": [float(jb), float(jbp)],
        "adf_residuos": [float(ar[0]), float(ar[1])],
        "robustez_lag": {"wald_F": sc(wl.statistic), "wald_p": sc(wl.pvalue),
                         "b_DxLIliq": float(ml.params["DxLI"]),
                         "p_DxLIliq": float(ml.pvalues["DxLI"])},
        "regimen_gradual": {"wald_F": sc(wg.statistic), "wald_p": sc(wg.pvalue),
                            "RxI": float(mg.params["RxI"]), "p_RxI": float(mg.pvalues["RxI"])},
    }
    json.dump(out, open("output/tablas/etapa5_revisada_tutor.json", "w"), indent=1)
    print("Guardado: etapa5_modelo_rp_centrado.csv + etapa5_revisada_tutor.json")


if __name__ == "__main__":
    main()
