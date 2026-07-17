# MANIFIESTO DE ENTREGA — versión canónica

**Fecha de congelamiento:** 17/07/2026 (rev. 3: corte principal 2011:T2, do-file y Word consistentes)
**Regla:** estos son los únicos archivos válidos para depósito, auditoría y defensa.
Cualquier otro archivo con nombre similar es una versión anterior archivada (carpeta `archivo_versiones/`).

| Archivo | SHA-256 |
|---|---|
| `Scavone_Fernandez_Tesis_FINAL_2026-07-17.docx` | `ab8d067047f2d0384c32ef13e4a1dbad5fea722b94aabf4453095508afab4f41` |
| `stata/tesis_cda.do` | `458a73b845a0dda56eaeb1c08e157e62ed51f9a015022fd316f80ce8a157cd08` |
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

## Rev. 3 — cerrada (migración al corte 2011T2; do-file y Word consistentes)
- **Decisión del autor (responde al hallazgo M3 del profesor):** el corte institucional
  principal pasa de 2011:T1 a **2011:T2** (trimestre real de la Resolución N.º 22, Acta N.º 31,
  del 18-may-2011). T1 queda como sensibilidad.
- **do-file (`stata/tesis_cda.do`) reescrito** sobre la estructura limpia revisada: rutas portables
  (`args root`), log, asserts de integridad (84 obs contiguas, sin missing en la serie principal,
  chequeo de fórmulas de inflación), corte 2011T2 principal, y ahora genera las **3 figuras nativas
  en Stata** (`Figura_8_1_DR`, `Figura_8_2_Fisher`, `Figura_8_3_Composicion`), más el quiebre
  estructural (ADF/PP/KPSS + Bai-Perron con `xtbreak`/`sbsingle`), la composición Shapley/LMG, la
  dominancia estandarizada y la cointegración de Engle-Granger.
- **HECHO — Word actualizado a T2 desde la corrida real en Stata 17** (log `output/tesis_cda.log`):
  β4 = +3,89 (p<0,001); Wald de pendientes F = 2,00 (p = 0,142); Chow F = 26,00 (p<0,001);
  Fisher Δrp = +3,89 y ΔInflación = −3,07 (ambas p<0,001); quiebre endógeno de la prima rp en
  2011:T2 (sup-Wald = 43,66; p<0,001). 92 cifras reemplazadas y verificadas: **ningún número T1
  quedó en el texto principal**; el 2011:T1 aparece solo como robustez de la fecha. do-file ↔ Word
  **consistentes** (condición C1 satisfecha).
- **Python retirado del set canónico:** la validación ya no depende de una réplica en Python; el
  do-file `stata/tesis_cda.do` se ejecutó en Stata 17 y produce cada cifra del documento (el log de
  esa corrida es la evidencia de reproducibilidad). El script Python anterior queda archivado.

## Rev. 2 (qué cambió respecto del congelamiento inicial)
- Números de inferencia HAC del Word alineados a la convención exacta de Stata (`F = 0,33; p = 0,717`
  para el Wald de pendientes; Chow `F = 15,55`; y análogos), coincidentes con la réplica independiente
  del informe de auditoría (F = 0,334; p = 0,717).
- Do-file: prueba t de Fisher sobre la misma muestra de la descomposición (n = 80); fechado
  Bai-Perron de la prima rp agregado (quiebre único: 2011:T2); comentarios actualizados.
- `base_cda.xlsx`: hojas `diccionario` y `fuentes_trazabilidad` reconstruidas (se habían perdido
  al integrarse la serie TPM); hoja de datos sin cambios.
