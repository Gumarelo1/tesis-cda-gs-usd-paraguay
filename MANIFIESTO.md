# MANIFIESTO DE ENTREGA — versión canónica

**Fecha de congelamiento:** 17/07/2026
**Regla:** estos son los únicos archivos válidos para depósito, auditoría y defensa.
Cualquier otro archivo con nombre similar es una versión anterior archivada.

| Archivo | SHA-256 |
|---|---|
| `Scavone_Fernandez_Tesis_FINAL_2026-07-17.docx` | `aa785d7bcc8cb677c927c00bbe79557ceb4a21651415eb044747dd87fab291f2` |
| `stata/tesis_cda.do` | `66fcd5d9e0302946d225f0c525d70741dc90fd54e0783a1895e37045f6cd5b43` |
| `data/processed/base_cda.csv` | `f5700221ca2ccab9389e0e3285b0cfc68d67864ee91bfd4b487ef4087b19631e` |
| `data/processed/base_cda.xlsx` | `423aaf6106e615c972ee72f9f78581bbc5bf8f93d014a3d94c92d1f9bdd97b3a` |
| `data/processed/base_cda.dta` | `5ffb7e9a7e73cbd82c7f6ddfe32e0961949ea1f1512493755fd904c54613284b` |
| `data/processed/sources.md` | `37812ebd9a4a47558e5bb8489054813b34c8a7be8d613af43eb1a51c2ee22237` |

## Correspondencia única (condición de la auditoría)
- **Serie principal:** tramo CDA ≤365 días, tasa efectiva, bancos (elección teórica: única fuente
  y definición homogénea 2004–2024, sin empalme). La ponderada queda solo como sensibilidad (Etapa 7).
- **Todas las tablas y cifras del Word** provienen de una única corrida del do-file sobre esta base
  (log en `output/tesis_cda.log` al ejecutar).
- **Datos crudos oficiales** en `data/raw/` (246 boletines de tasas, Anexo Estadístico, TPM.xlsx, FRED).
- **Trazabilidad:** hoja `diccionario` y `fuentes_trazabilidad` dentro de `base_cda.xlsx`, y `data/processed/sources.md`.
