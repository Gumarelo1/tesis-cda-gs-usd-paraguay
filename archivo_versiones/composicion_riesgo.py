# -*- coding: utf-8 -*-
"""
Composicion interna de la prima de riesgo (rp) — equivalente Python de la
ETAPA 5-BIS del do-file (stata/tesis_cda.do). Genera:
  - output/tablas/etapa4b_composicion_riesgo.csv          (estatica pre/post)
  - output/tablas/etapa4c_composicion_riesgo_ampliada.csv (Shapley con globales)
  - output/tablas/etapa4d_evolucion_composicion_riesgo*.csv (ventanas 12 y 20T)
  - output/figuras/etapa4d_evolucion_composicion_riesgo.png
HONESTIDAD: analisis descriptivo/exploratorio; se reporta tal cual sale.
"""
import math
from itertools import combinations

import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "data/processed/base_cda.csv"


def cargar() -> pd.DataFrame:
    df = pd.read_csv(BASE).sort_values(["anio", "trimestre"]).reset_index(drop=True)
    df["fecha_dt"] = pd.PeriodIndex(
        df["anio"].astype(str) + "Q" + df["trimestre"].astype(str), freq="Q"
    ).to_timestamp()
    df["difn"] = df["i_Gs"] - df["i_USD"]
    df["dep"] = df["TC"].pct_change() * 100
    df["infl_py"] = (df["IPC_PY"] / df["IPC_PY"].shift(4) - 1) * 100
    df["infl_us"] = (df["IPC_US"] / df["IPC_US"].shift(4) - 1) * 100
    df["InflDiff"] = df["infl_py"] - df["infl_us"]
    df["rp"] = df["difn"] - df["InflDiff"]
    df["VolTC"] = df["dep"].rolling(8, min_periods=2).std()
    df["Iliq"] = df["dolariz"]
    df["D2011"] = (df["anio"] >= 2011).astype(int)
    df["FF"] = df["FEDFUNDS"]
    df["dFF"] = df["FEDFUNDS"].diff()
    df["Crisis0809"] = (
        ((df["anio"] == 2008) & (df["trimestre"] >= 2))
        | ((df["anio"] == 2009) & (df["trimestre"] <= 3))
    ).astype(int)
    return df


def lmg2(y, x1, x2):
    """Shapley/LMG exacto para 2 regresores (formula cerrada)."""
    r1 = sm.OLS(y, sm.add_constant(np.column_stack([x1]))).fit().rsquared
    r2 = sm.OLS(y, sm.add_constant(np.column_stack([x2]))).fit().rsquared
    r12 = sm.OLS(y, sm.add_constant(np.column_stack([x1, x2]))).fit().rsquared
    return 0.5 * (r1 + r12 - r2), 0.5 * (r2 + r12 - r1), r12


def lmg_k(y, X: pd.DataFrame):
    """Shapley/LMG exacto para k regresores (2^k subconjuntos)."""
    names = list(X.columns)
    cache: dict = {}

    def r2_of(s):
        if not s:
            return 0.0
        k = tuple(sorted(s))
        if k not in cache:
            cache[k] = sm.OLS(y, sm.add_constant(X[list(k)])).fit().rsquared
        return cache[k]

    n = len(names)
    contrib = {}
    for v in names:
        others = [o for o in names if o != v]
        tot = 0.0
        for r in range(n):
            for sub in combinations(others, r):
                w = math.factorial(r) * math.factorial(n - r - 1) / math.factorial(n)
                tot += w * (r2_of(set(sub) | {v}) - r2_of(set(sub)))
        contrib[v] = tot
    return contrib, r2_of(set(names))


def main() -> None:
    df = cargar()
    d = df.dropna(subset=["rp", "VolTC", "Iliq"]).reset_index(drop=True)

    # (a) estatica pre/post
    rows = []
    for lbl, sub in [("PRE-2011", d[d.D2011 == 0]), ("POST-2011", d[d.D2011 == 1]),
                     ("TOTAL 2004-2024", d)]:
        c1, c2, r2 = lmg2(sub["rp"], sub["VolTC"], sub["Iliq"])
        tot = c1 + c2
        rows.append({"label": lbl, "n": len(sub), "r2_vol": None, "r2_iliq": None,
                     "r2_both": r2, "corr_vol": sub["rp"].corr(sub["VolTC"]),
                     "corr_iliq": sub["rp"].corr(sub["Iliq"]),
                     "pct_vol": 100 * c1 / tot, "pct_iliq": 100 * c2 / tot})
        print(f"(a) {lbl}: n={len(sub)} R2={r2:.3f} "
              f"%VolTC={100*c1/tot:.1f} %Iliq={100*c2/tot:.1f}")
    pd.DataFrame(rows).to_csv("output/tablas/etapa4b_composicion_riesgo.csv", index=False)

    # (b) ampliada con factores globales (Shapley exacto, periodo completo)
    cands = ["VolTC", "Iliq", "FF", "dFF", "Crisis0809", "dIliq"]
    df["dIliq"] = df["dolariz"].diff()
    db = df.dropna(subset=["rp"] + cands)
    contrib, r2 = lmg_k(db["rp"], db[cands])
    amp = [{"periodo": "TOTAL", "variable": v, "contrib_R2": c,
            "pct_del_explicado": 100 * c / r2,
            "corr_con_rp": db["rp"].corr(db[v])} for v, c in contrib.items()]
    amp.append({"periodo": "TOTAL", "variable": "__R2_TOTAL__", "contrib_R2": r2,
                "pct_del_explicado": 100.0, "corr_con_rp": np.nan})
    pd.DataFrame(amp).to_csv(
        "output/tablas/etapa4c_composicion_riesgo_ampliada.csv", index=False)
    print(f"(b) R2 total={r2:.3f}; " + "; ".join(
        f"{v}={100*c/r2:.1f}%" for v, c in contrib.items()))

    # (c) evolucion en ventanas moviles 12T y 20T + figura
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    for ax, win in zip(axes, [12, 20]):
        rws = []
        for i in range(win - 1, len(d)):
            w = d.iloc[i - win + 1:i + 1]
            if w["rp"].std() < 1e-8:
                continue
            c1, c2, r2w = lmg2(w["rp"].values, w["VolTC"].values, w["Iliq"].values)
            tot = c1 + c2
            rws.append({"fecha": w["fecha_dt"].iloc[-1],
                        "pct_VolTC": 100 * c1 / tot if tot > 1e-9 else np.nan,
                        "pct_Iliq": 100 * c2 / tot if tot > 1e-9 else np.nan,
                        "R2": r2w})
        e = pd.DataFrame(rws)
        suf = "" if win == 12 else "_win20"
        e.to_csv(f"output/tablas/etapa4d_evolucion_composicion_riesgo{suf}.csv",
                 index=False)
        ax.stackplot(e.fecha, e.pct_VolTC, e.pct_Iliq,
                     labels=["Riesgo cambiario (VolTC)", "Liquidez/dolarizacion (Iliq)"],
                     colors=["#c0392b", "#2980b9"], alpha=0.8)
        ax.axvline(pd.Timestamp("2011-01-01"), color="black", lw=2, ls="--")
        ax.axvspan(pd.Timestamp("2008-04-01"), pd.Timestamp("2009-09-01"),
                   color="gray", alpha=0.25)
        ax.axvspan(pd.Timestamp("2020-01-01"), pd.Timestamp("2020-12-01"),
                   color="gray", alpha=0.25)
        ax.set_ylim(0, 100)
        ax.set_ylabel(f"% de la prima de riesgo\n(ventana {win}T)")
        ax.set_title(f"Composicion de la prima de riesgo — ventana movil {win} trimestres")
    axes[0].legend(loc="upper left", fontsize=9)
    axes[1].set_xlabel("Trimestre (linea negra=2011 | gris=crisis 2008-09 y pandemia 2020)")
    plt.tight_layout()
    plt.savefig("output/figuras/etapa4d_evolucion_composicion_riesgo.png", dpi=140)
    print("(c) figura y tablas de evolucion guardadas")

    # (d) CONTRASTE FORMAL de la hipotesis (forma general):
    #   H0: la composicion de la prima de riesgo NO cambio tras 2011.
    #   Se contrasta sobre rp (objeto de la teoria), con interacciones de
    #   regimen y errores HAC-4. Equivalente al bloque (d) del do-file.
    dd = d.copy()
    dd["D_VolTC"] = dd["D2011"] * dd["VolTC"]
    dd["D_Iliq"] = dd["D2011"] * dd["Iliq"]
    X = sm.add_constant(dd[["VolTC", "Iliq", "D2011", "D_VolTC", "D_Iliq"]])
    m = sm.OLS(dd["rp"], X).fit(cov_type="HAC", cov_kwds={"maxlags": 4}, use_t=True)
    w_slopes = m.wald_test("D_VolTC = 0, D_Iliq = 0", use_f=True, scalar=True)
    w_full = m.wald_test("D2011 = 0, D_VolTC = 0, D_Iliq = 0", use_f=True, scalar=True)
    filas = [{"parametro": v, "coef": m.params[v], "SE_HAC": m.bse[v],
              "p": m.pvalues[v]} for v in ["VolTC", "Iliq", "D2011", "D_VolTC", "D_Iliq"]]
    filas.append({"parametro": "F_test_pendientes(D_VolTC=D_Iliq=0)",
                  "coef": float(w_slopes.statistic), "SE_HAC": np.nan,
                  "p": float(w_slopes.pvalue)})
    filas.append({"parametro": "F_test_Chow(D2011=D_VolTC=D_Iliq=0)",
                  "coef": float(w_full.statistic), "SE_HAC": np.nan,
                  "p": float(w_full.pvalue)})
    pd.DataFrame(filas).to_csv(
        "output/tablas/etapa4e_test_cambio_composicion.csv", index=False)
    print(f"(d) test formal cambio de composicion: "
          f"F={float(w_slopes.statistic):.3f} p={float(w_slopes.pvalue):.4f} "
          f"(Chow completo: F={float(w_full.statistic):.3f} "
          f"p={float(w_full.pvalue):.4f}); D_Iliq={m.params['D_Iliq']:+.4f} "
          f"p={m.pvalues['D_Iliq']:.4f}")


if __name__ == "__main__":
    main()
