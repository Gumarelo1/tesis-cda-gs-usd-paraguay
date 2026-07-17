# MANIFIESTO DE ENTREGA — versión canónica

**Fecha de congelamiento:** 17/07/2026 (rev. 2: números de inferencia alineados a Stata)
**Regla:** estos son los únicos archivos válidos para depósito, auditoría y defensa.
Cualquier otro archivo con nombre similar es una versión anterior archivada (carpeta `archivo_versiones/`).

| Archivo | SHA-256 |
|---|---|
| `Scavone_Fernandez_Tesis_FINAL_2026-07-17.docx` | `ca11e6554c9d909f5c71632eeb0df79933e43c9de7e85cd3dafc91afb2851560` |
| `stata/tesis_cda.do` | `9a98a77226af0fd01e2919bfeabcfaf9439ddf53a3017024f838e0823ce649fa` |
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

## Rev. 2 (qué cambió respecto del congelamiento inicial)
- Números de inferencia HAC del Word alineados a la convención exacta de Stata (`F = 0,33; p = 0,717`
  para el Wald de pendientes; Chow `F = 15,55`; y análogos), coincidentes con la réplica independiente
  del informe de auditoría (F = 0,334; p = 0,717).
- Do-file: prueba t de Fisher sobre la misma muestra de la descomposición (n = 80); fechado
  Bai-Perron de la prima rp agregado (quiebre único: 2011:T2); comentarios actualizados.
- `base_cda.xlsx`: hojas `diccionario` y `fuentes_trazabilidad` reconstruidas (se habían perdido
  al integrarse la serie TPM); hoja de datos sin cambios.
