# Fuentes y trazabilidad — Base trimestral CDA Gs/USD (2004Q1–2024Q4)

> ⚠️ **ESTADO VIGENTE DE LA BASE (leer primero).** La serie DEFINITIVA de tasas es el
> **tramo "Certif. Dep. de Ahorro ≤365 días" (tasa efectiva, Bancos)** de «Promedio de
> Tasas» (SB/BCP). Con esa serie **la base está COMPLETA: los 84 trimestres tienen valor
> de tasa CDA en ambas monedas, sin ningún hueco trimestral** (verificado: 0 faltantes en
> `i_Gs` e `i_USD`). Detalle de construcción y meses en §9. Las **secciones 2 y 3 de abajo
> describen la serie PONDERADA ANTERIOR** (`i_*_pond`, conservada solo para sensibilidad),
> que sí tenía celdas vacías en 2006Q4/2010Q4; **NO aplican a la serie definitiva ≤365** y
> se conservan únicamente por trazabilidad del proceso.

**Tesis:** "Diferencial de Retornos entre Certificados de Depósito en Guaraníes y Dólares en el Sistema Bancario Paraguayo, 2004–2024" (Laviero Scavone, Matías Fernández — UC Asunción).
**Fecha de extracción:** 2026-06-30.
**Ventana:** 2004Q1–2024Q4 = 84 trimestres.
**Regla de oro:** ningún dato inventado, estimado ni interpolado. Todo valor se rastrea a un archivo oficial descargado en `data/raw/`. Los huecos se dejan vacíos y se declaran abajo.

---

## 1. Archivos crudos descargados

| Archivo local | Fuente | URL exacta |
|---|---|---|
| `data/raw/bcp/Anexo_Estadistico_Informe_Economico.xlsx` | BCP – Anexo Estadístico del Informe Económico (planilla histórica, 94 hojas) | https://www.bcp.gov.py/documents/20117/213069/Anexo_Estad%C3%ADstico_del_Informe_Econ%C3%B3mico_30_06_2026.xlsx/b62faece-95c4-5cbc-429c-01696c19c4f3?t=1782850649323 |
| `data/raw/fred/FEDFUNDS.csv` | FRED (Federal Reserve St. Louis) | https://fred.stlouisfed.org/graph/fredgraph.csv?id=FEDFUNDS |
| `data/raw/fred/CPIAUCSL.csv` | FRED | https://fred.stlouisfed.org/graph/fredgraph.csv?id=CPIAUCSL |
| `data/raw/bcp/if_pdfs/*.pdf` (16 PDFs) | BCP – Informes de Indicadores Financieros mensuales (2004–2012) | ver `data/raw/bcp/if_pdfs/_INVENTARIO.csv` |

También se conservan páginas HTML/renderizados de reconocimiento en `data/raw/bcp/pages/` y `data/raw/bcp/rendered/` (no son fuentes de dato, solo evidencia del proceso de descubrimiento de URLs).

---

## 2. Detalle por variable

### `i_Gs` — Tasa efectiva pasiva CDA en guaraníes (%)
- **2011-01 → 2024:** Anexo, hoja **CUADRO 31** (Moneda Nacional), columna **E "CDA"**. Fecha en col A.
- **2004-2010:** Informes de Indicadores Financieros mensuales (PDF). En cada PDF, página **"Tasas bancarias efectivas (promedios mensuales)"**, sección **"Pasivas en M/N"**, fila **"- C D A"** (o "- Cert. Depósitos Ahorro" en formato viejo). Es la tasa del **sistema bancario** (no financieras). Traza mes-a-mes en `data/processed/cda_if_pdfs_crudo.csv`.
- **Agregación:** trimestre = **promedio de los 3 meses**. Variante `i_Gs_fdt` = valor del **último mes** del trimestre.
- **Cobertura:** 81/84 (promedio) — faltan 2005Q2, 2006Q4, 2010Q4 (ver §3). `i_Gs_fdt`: 83/84 (falta solo 2010Q4).

### `i_USD` — Tasa efectiva pasiva CDA en dólares (%)
- **2011-01 → 2024:** Anexo, hoja **CUADRO 31 (Cont.)** (Moneda Extranjera), columna **E "CDA"**.
- **2004-2010:** ídem `i_Gs` pero sección **"Pasivas en M/E"**, fila "- C D A".
- **Agregación:** promedio de 3 meses (`i_USD_fdt` = último mes).
- **Cobertura:** 82/84 — faltan 2006Q4, 2010Q4. `i_USD_fdt`: 83/84.

### `TC` — Tipo de cambio nominal Gs/USD
- Anexo, hoja **CUADRO 60a** "Tipo de cambio nominal del guaraní", fecha en **col B**, valor en **col E "USD"**.
- **IMPORTANTE:** la fuente publica **promedio mensual**, no cotización de fin de día. La regla del proyecto pide "último día hábil del trimestre"; esa serie diaria de cierre **no está disponible en fuente oficial histórica** para 2004-2011 (la webapp de cotización referencial fluctuante del BCP solo llega a 2012 y usa otra metodología — mercado fluctuante con clientes no financieros). Por consistencia en toda la ventana se usa el **promedio mensual del mes de cierre** del trimestre. **Decisión para el auditor.**
- **Agregación:** valor del mes de cierre. **Cobertura:** 84/84.

### `dolariz` — Ratio de dolarización de depósitos (%)
- Anexo, hoja **CUADRO 23a** "Depósitos del sector privado en bancos y financieras", **col O "Part. %"** de moneda extranjera (fracción × 100). Cuadro 23b (cooperativas) NO usado.
- **Agregación:** valor del mes de cierre. **Cobertura:** 84/84. (2024Q2–Q4 provienen de filas marcadas "* Cifras preliminares" en la fuente.)

### `IPC_PY` — IPC Paraguay, nivel (base dic-2017=100, serie empalmada)
- Anexo, hoja **CUADRO 15**, **col V "Índice general"**. **Agregación:** nivel del mes de cierre. **Cobertura:** 84/84.

### `infl_PY_yoy` — Inflación interanual Paraguay (%)
- Anexo, **CUADRO 15**, **col Y "Inflación total – Interanual"** (publicada por el BCP; NO recalculada). Mes de cierre. **Cobertura:** 84/84.

### `IPC_US` — CPIAUCSL, nivel
- FRED, serie CPIAUCSL. Nivel del mes de cierre. **Cobertura:** 84/84.

### `infl_US_yoy` — Inflación interanual EEUU (%)
- Derivada de CPIAUCSL: `(CPI_t / CPI_{t-12} − 1) × 100`. Transformación determinista del nivel oficial (no es dato inventado). Mes de cierre. **Cobertura:** 84/84.

### `FEDFUNDS` — Federal Funds Effective Rate (%)
- FRED, serie FEDFUNDS. **Agregación:** promedio de 3 meses (consistente con tasas). **Cobertura:** 84/84.

### `i_Gs_vista`, `i_USD_vista` — Tasas "a la vista" por moneda (REFERENCIA)
- Anexo, CUADRO 31 / 31 (Cont.), **col C "A la vista"**. NO son CDA; se incluyen solo como referencia auxiliar. Promedio de 3 meses. Cobertura 84/84.

---

## 3. Huecos de la serie PONDERADA ANTERIOR (superseded — NO aplican a la ≤365)

> Esta sección describe la serie ponderada vieja (`i_*_pond`). **La serie definitiva ≤365
> no tiene estos huecos: está completa 84/84** (ver banner arriba y §9).


| Trimestre | Variable(s) | Motivo |
|---|---|---|
| **2005Q2** | `i_Gs` | El mes **abr-2005** es **dudoso** (dos IF PDFs de 2005 leen 13.48%; el IF de dic-2005 lee 8.89% — aparente revisión del BCP). Al no haber consenso, abr-2005 se deja NaN → el promedio trimestral no se computa. Detalle en `cda_conflictos.csv`. `i_Gs_fdt` (jun-2005) sí tiene dato. |
| **2006Q4** | `i_Gs`, `i_USD` | **nov-2006** no extraíble: el PDF (`nov (2).pdf`) es un **escaneo sin capa de texto**. No se hace OCR para no arriesgar fabricación. (dic-2006 sí está, vía IF dic-2007 → `*_fdt` de 2006Q4 tiene dato.) |
| **2010Q4** | `i_Gs`, `i_USD`, `*_fdt` | oct/nov/dic-2010 no localizados: el Anexo (Cuadro 31) recién arranca CDA en 2011-01, y no se halló un IF mensual público que cubra oct–dic 2010 (el último IF viejo accesible llega a sep-2010; el repositorio DSpace `repositorio.bcp.gov.py` está inaccesible desde el entorno). |

**Ningún otro hueco.** El resto de variables tiene 84/84.

---

## 4. Nota metodológica sobre la serie CDA (empalme 2010/2011)

`i_Gs`/`i_USD` combinan dos publicaciones del BCP, ambas "sistema bancario, CDA, tasa efectiva":
- **2004-2010:** Informes de Indicadores Financieros mensuales (tabla "Tasas bancarias efectivas").
- **2011-2024:** Anexo Estadístico, Cuadro 31.

En el empalme hay un **salto real** (Gs CDA ≈ 7.7% en 2010Q3 → ≈ 11.1% en 2011Q1) que corresponde al **endurecimiento monetario efectivo del BCP en 2011** (suba de tasas), no a un error de extracción ni a un cambio de definición conocido. Aun así, al provenir de dos publicaciones distintas, **se advierte al auditor/analista**: conviene verificar si existe una discontinuidad de definición. Los valores de cada tramo son fieles a su fuente.

Verificación cruzada 2004-2010: de ~160 pares mes-moneda con solapamiento entre PDFs, **solo 1 discrepancia** > 0.15 p.p. (abr-2005, arriba). El resto coincide.

---

## 5. Archivos intermedios (para que el analista pueda re-agregar)

- `data/processed/i_Gs_mensual.csv`, `i_USD_mensual.csv` — CDA y "a la vista" mensuales (anexo).
- `data/processed/IPC_PY_mensual.csv`, `dolariz_mensual.csv`, `TC_prom_mensual.csv` — anexo.
- `data/processed/cda_if_pdfs_crudo.csv` — **traza fila-por-fila** de la CDA leída de cada PDF (archivo, moneda, año, mes, valor).
- `data/processed/cda_conflictos.csv` — meses con lecturas discrepantes entre PDFs (dudosos).
- `data/raw/bcp/if_pdfs/_INVENTARIO.csv` — mapeo de cada PDF a su ventana de contenido y URL.
- `data/processed/LIBRO_DE_FUENTES.csv` — trazabilidad estructurada por variable.

---

## 6. Fuentes NO usadas como dato primario (por regla del proyecto)

- **Banco Mundial / FMI / OCDE**: prohibidas como primaria. (FRED `PRYCPIALLMINMEI` se descartó: además el endpoint devolvió HTML de error).
- **repositorio.bcp.gov.py** (DSpace con IF históricos 2004-2010 sistemáticos): **inaccesible** desde el entorno (timeout/conexión rechazada, tanto por curl como por navegador headless). Habría sido la vía ideal para cerrar los 4 meses faltantes; se deja constancia.

---

## 7. Actualización post-auditoría (nota de cobertura)

**2005Q2 `i_Gs` = 9.03 (recuperado).** La auditoría resolvió el conflicto de lecturas de abr-2005 (8.89 adoptado; ver `AUDITORIA.md` §4a y `CORRECCIONES_auditor.csv`). En la **serie ponderada anterior** la cobertura de `i_Gs` era 82/84 (huecos en 2006Q4 y 2010Q4). **En la serie definitiva ≤365 la cobertura es 84/84, sin huecos** (ver banner al inicio y §9). Las secciones 2-3 reflejan el estado PRE-auditoría de la serie ponderada y se conservan por trazabilidad.

## 8. RECUPERACIÓN DEL HUECO 2010Q4 (2026-07-13)

Los IF mensuales de **oct/nov/dic-2010** fueron recuperados del BCP vía URL directa
(no estaban listados en ninguna página, pero existen en el repositorio institucional):
- https://www.bcp.gov.py/documents/d/institucional/if_octubre_2010-pdf
- https://www.bcp.gov.py/documents/d/institucional/if_noviembre_2010-pdf
- https://www.bcp.gov.py/documents/d/institucional/if_diciembre_2010-pdf
(archivos en `data/raw/bcp/if_pdfs/hunt/`; tabla "Tasa de Interés Pasiva" MN pág. 4 y ME pág. 5, fila CDA)

Valores CDA (misma definición que el resto de la serie 2004-2010; validado: sep-10 de estos
PDFs = 7,68 Gs / 4,44 USD, idéntico a la serie ya auditada): Gs 7,73/7,78/8,79; USD 4,63/4,50/4,65.
→ 2010Q4: i_Gs=8,10 (promedio), i_USD=4,5933, fdt=8,79/4,65.
Verificación adicional: los 6 valores caen dentro de los rangos por tramos de plazo publicados
en "Promedio de Tasas Mensuales" (ver docs/08). **Hueco restante: solo 2006Q4** (nov-2006
irrecuperable: PDF escaneado sin texto en un canal, recorte de prensa en otro, sin URL directa).

## 9. CAMBIO DE SERIE PRINCIPAL DE TASAS (2026-07-13, decisión del autor)

Por decisión del autor, la serie principal de tasas CDA pasa a ser el tramo
**"CERTIFIC.DEP. DE AHORRO <= 365 DIAS" (bancos, tasa EFECTIVA promedio ponderado)**,
de los boletines mensuales "Promedio de Tasas" de la Superintendencia de Bancos:
https://www.bcp.gov.py/promedio-de-tasas-mensuales
- M/L = guaraníes (i_Gs); M/E = dólares (i_USD). 242 boletines mensuales descargados y
  parseados (`data/raw/bcp/promedio_tasas/full/`, crudo mes a mes en
  `data/processed/cda365_mensual_crudo.csv`, trimestral en `cda365_trimestral.csv`).
- VENTAJA: una sola fuente y una sola definición para TODO 2004-2024 (elimina el empalme
  IF/Cuadro 31 y cubre los 84 trimestres, incluidos 2006Q4 y 2010Q4).
- 10 meses NO publicados en esa sección (jul-2005, mar-2006, nov-2006, dic-2006, ene-2007,
  may-2007, jul-2007, abr-2013, may-2013, abr-2014; cuatro de ellos son recortes de prensa
  subidos por error por el BCP). Los 8 trimestres afectados promedian los meses disponibles
  (columna `n_meses_365` de la base); i_*_fdt usa el último mes disponible del trimestre.
  NO se imputó ningún valor.
- La serie consolidada anterior (IF 2004-2010 + Cuadro 31 2011-2024) se conserva en la base
  como `i_Gs_cons`/`i_USD_cons`(+`_fdt`) para robustez. Correlación entre ambas: 0,86 (Gs) /
  0,86 (USD). El "salto" 2010Q4→2011Q1 de la serie consolidada (8,10→11,10, en parte
  definicional por el empalme) no existe en la serie homogénea (6,79→6,88): el aumento real
  ocurre gradualmente hacia 2011Q2 (9,45).

## 9. CAMBIO DE SERIE PRIMARIA CDA (2026-07-13, pedido por el autor)

La serie primaria `i_Gs`/`i_USD` pasa a ser el tramo **"CERTIFIC. DEP. DE AHORRO <= 365 DIAS"
(tasa EFECTIVA, BANCOS, M/N y M/E)** de los boletines mensuales "Promedio de Tasas"
de la Superintendencia de Bancos (https://www.bcp.gov.py/promedio-de-tasas-mensuales),
2004-2024 completo. Ventajas: (i) UNA sola fuente y UNA sola definición para toda la
muestra -> desaparece el empalme 2010/2011; (ii) tramo de plazo fijo -> sin sesgo de
composición por cambios en la estructura de plazos.

- Crudos: `data/raw/bcp/promedio_tasas/full/` (246 PDF) + `cda365_mensual_crudo.csv` (242 meses).
- Meses NO publicados por el BCP (verificado contra su propio listado + sondas de URL directa):
  2005-07, 2006-03, 2006-11*, 2006-12*, 2007-01*, 2007-05*, 2007-07, 2013-04, 2013-05, 2014-04
  (*el archivo publicado es un recorte de prensa sin la tabla).
- Convención trimestral: promedio de los meses publicados; `meses_cobertura` y
  `ultimo_mes_publicado` en la base documentan cada trimestre. 8 trimestres con
  cobertura parcial (2006Q4 y 2013Q2 con 1 mes; declarados, sin imputar).
- `i_*_fdt` = último mes publicado del trimestre (en 2006Q1 y 2006Q4 no es el mes de cierre).
- La serie anterior (promedio ponderado consolidado IF 2004-2010 + Cuadro 31 2011-2024)
  se conserva como `i_*_pond` para contraste; corr mensual nuevo-viejo: Gs 0,87 / USD 0,83.
