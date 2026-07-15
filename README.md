# Recomposición de la Prima de Riesgo en el Diferencial de Retornos CDA Gs/USD

**Sistema Bancario Paraguayo, 2004–2024**

Tesis de grado — Laviero Scavone y Matías Fernández (Universidad Católica de Asunción).
Este repositorio contiene la **base de datos** y el **do-file de Stata** para auditoría.

---

## Contenido

| Archivo | Descripción |
|---|---|
| `Scavone_Fernandez_Entrega_FINAL_365.docx` | Documento de la tesis (versión final). |
| `stata/tesis_cda.do` | Do-file de Stata 17. Implementa las 7 etapas del análisis econométrico. |
| `data/processed/base_cda.csv` | Base trimestral, 84 observaciones (2004T1–2024T4). Formato principal. |
| `data/processed/base_cda.xlsx` | Misma base, formato Excel. |
| `data/processed/base_cda.dta` | Misma base, formato nativo de Stata. |
| `data/processed/sources.md` | Trazabilidad completa: de qué archivo oficial y qué celda sale cada dato. |

Los tres formatos de la base son idénticos en contenido.

## La base de datos

- **84 trimestres completos** (2004T1–2024T4), sin huecos trimestrales.
- **Serie de tasas CDA:** tramo "Certif. Dep. de Ahorro ≤ 365 días" (tasa efectiva, bancos)
  de los boletines mensuales *"Promedio de Tasas"* de la Superintendencia de Bancos del BCP.
  Fuente única y una sola definición para todo el período (sin empalme).
- **Otras fuentes oficiales:** tipo de cambio, dolarización de depósitos e IPC de Paraguay
  (Anexo Estadístico del Informe Económico, BCP); IPC de EE.UU. y tasa de fondos federales (FRED).
- **Regla de datos:** ningún valor inventado, estimado ni interpolado. Todo dato se rastrea a
  un archivo oficial. Ver `data/processed/sources.md` para la trazabilidad celda por celda.

## El do-file (`stata/tesis_cda.do`)

Estructura (7 etapas):

1. **Etapa 1** — series de retorno y estadística descriptiva.
2. **Etapa 2** — coeficiente de Fama (paridad descubierta) y prueba conservadora.
3. **Etapa 3** — estacionariedad (ADF/KPSS/Phillips-Perron) y quiebre estructural (Bai-Perron).
4. **Etapa 4** — descomposición de Fisher (prima de inflación vs. prima de riesgo).
5. **Etapa 5** — modelo MCO-HAC con interacciones sobre la prima de riesgo (núcleo), más la
   composición interna (Shapley/LMG) y la revisión metodológica (centrado, dominancia,
   integración, Jarque-Bera, cointegración de Engle-Granger, régimen gradual).
6. **Etapa 6** — diagnósticos y robustez (incluye control por la crisis de 2008-09 y
   sensibilidad al rezago HAC).
7. **Etapa 7** — sensibilidad a la definición de la serie de tasas.

Nota metodológica: el modelo se estima sobre la **prima de riesgo** `rp = diferencial de tasas −
diferencial de inflación`, con la inflación **restada** (nunca como regresor), para evitar la
circularidad de Fisher.

### Cómo correrlo

1. Abrir `stata/tesis_cda.do` y ajustar **una sola línea** — la ruta de la carpeta:
   ```stata
   global root "C:/ruta/a/esta/carpeta"
   ```
2. Requiere Stata 17 y los módulos comunitarios `rangestat`, `estout`, `xtbreak` y `kpss`
   (el do-file los instala automáticamente desde SSC si faltan y hay conexión).
3. Ejecutar: `do "stata/tesis_cda.do"`.

## Reproducibilidad

Todos los coeficientes, errores estándar y p-valores se reportan tal como salen; no se ajustó
ningún dato ni especificación para forzar la hipótesis.
