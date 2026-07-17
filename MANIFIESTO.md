# MANIFIESTO DE ENTREGA — versión canónica

**Fecha de congelamiento:** 17/07/2026 (rev. 2: números de inferencia alineados a Stata)
**Regla:** estos son los únicos archivos válidos para depósito, auditoría y defensa.
Cualquier otro archivo con nombre similar es una versión anterior archivada (carpeta `archivo_versiones/`).

| Archivo | SHA-256 |
|---|---|
| `Scavone_Fernandez_Tesis_FINAL_2026-07-17.docx` | `2f2c1b25850f684e56388f91874d6f91bc7dcfbac0fe283377011bc78fc9e78e` |
| `stata/tesis_cda.do` | (rev.3 EN CURSO — ver nota abajo) |
| `python/corrida_canonica.py` | `95ec53c4dec7db70a7ef812cc90f3ccecf6f78f881a4e95cdfb8fc170fcd8012` |
| `data/processed/base_cda.csv` | `f5700221ca2ccab9389e0e3285b0cfc68d67864ee91bfd4b487ef4087b19631e` |
| `data/processed/base_cda.xlsx` | `7f2d7d0096d75fa9b9dae65736970986f7d2ac26bccb6b05a5d0e386e7c48f75` |
| `data/processed/base_cda.dta` | `5ffb7e9a7e73cbd82c7f6ddfe32e0961949ea1f1512493755fd904c54613284b` |
| `data/processed/sources.md` | `37812ebd9a4a47558e5bb8489054813b34c8a7be8d613af43eb1a51c2ee22237` |

## Correspondencia única (condición de la auditoría)
- **Serie principal:** tramo CDA ≤365 días, tasa efectiva, bancos (elección teórica: única fuente
  y definición homogénea 2004–2024, sin empalme). La ponderada queda solo como sensibilidad (Etapa 7).
- **Todas las tablas y cifras del Word** provienen de una única corrida del do-file sobre esta base
  (log en `output/tesis_cda.log` al ejecutar). La inferencia del documento usa la convención de
  Stata `newey` (ajuste de muestra pequeña y distribución t), de modo que la corrida en Stata 17
  reproduce las cifras del texto; la validación cruzada independiente (`python/corrida_canonica.py`,
  con la misma convención) las reproduce en `output/canonica/resultados.json`.
- **Datos crudos oficiales** en `data/raw/` (boletines de tasas, Anexo Estadístico, TPM.xlsx, FRED).
- **Trazabilidad:** hojas `diccionario` y `fuentes_trazabilidad` dentro de `base_cda.xlsx`,
  y `data/processed/sources.md`.

## Rev. 3 — EN CURSO (migración al corte 2011T2)
- **Decisión del autor (responde al hallazgo M3 del profesor):** el corte institucional
  principal pasa de 2011:T1 a **2011:T2** (trimestre real de la Resolución N.º 22, Acta N.º 31,
  del 18-may-2011). T1 queda como sensibilidad.
- **do-file (`stata/tesis_cda.do`) reescrito** sobre la estructura limpia revisada: rutas portables
  (`args root`), log, asserts de integridad (84 obs contiguas, sin missing en la serie principal,
  chequeo de fórmulas de inflación), corte 2011T2 principal, y ahora genera las **3 figuras nativas
  en Stata** (`Figura_8_1_DR`, `Figura_8_2_Fisher`, `Figura_8_3_Composicion`), más el quiebre
  estructural (ADF/PP/KPSS + Bai-Perron con `xtbreak`/`sbsingle`), la composición Shapley/LMG, la
  dominancia estandarizada y la cointegración de Engle-Granger.
- **PENDIENTE (para cerrar la rev.3 y restaurar la consistencia C1):** actualizar los números del
  Word de T1 a T2 (β4 = +3,89; Wald F = 2,00 p = 0,142; Chow F = 26,0; Fisher Δrp = +3,89, ΔInfl =
  −3,07; y análogos). Se hará **a partir de la corrida del do-file en Stata 17** (log
  `output/tesis_cda.log`), para que cada cifra del texto sea exactamente la que produce Stata.
  Hasta entonces, el Word conserva las cifras T1 y el do-file usa T2: es un estado de trabajo
  declarado, no una versión de entrega.

## Rev. 2 (qué cambió respecto del congelamiento inicial)
- Números de inferencia HAC del Word alineados a la convención exacta de Stata (`F = 0,33; p = 0,717`
  para el Wald de pendientes; Chow `F = 15,55`; y análogos), coincidentes con la réplica independiente
  del informe de auditoría (F = 0,334; p = 0,717).
- Do-file: prueba t de Fisher sobre la misma muestra de la descomposición (n = 80); fechado
  Bai-Perron de la prima rp agregado (quiebre único: 2011:T2); comentarios actualizados.
- `base_cda.xlsx`: hojas `diccionario` y `fuentes_trazabilidad` reconstruidas (se habían perdido
  al integrarse la serie TPM); hoja de datos sin cambios.
