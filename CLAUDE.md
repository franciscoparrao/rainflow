# rainflow — Modelos hidrológicos conceptuales en Rust ("airGR/HBV moderno")

> **Estado:** EN DESARROLLO (v0.1 iniciada 2026-06-10). GR4J + métricas implementados,
> paridad numérica con airGR verificada (max diff 6e-7 mm, test de regresión en
> `crates/rainflow-core/tests/airgr_parity.rs`).
> Familia de motores Rust del autor: SurtGIS, Hydroflux, Smelt, Anvil, Cantus, Criterium.
> Doc madre: `~/proyectos/ideas-motores-rust.md` (idea B1).

## Qué es
Motor de modelos lluvia-escorrentía agregados/semi-distribuidos con calibración
automática y métricas de bondad de ajuste. Operacional y rápido.

## El gap que llena
**Hydroflux** es un solver físico-distribuido (caro, alta resolución). Falta el
otro extremo: modelos **conceptuales rápidos** (GR4J, HBV, Sacramento) para
pronóstico operacional y corridas masivas. Hoy: airGR/HBV (R), TUWmodel.

## Alcance MVP (v0.1)
- [x] GR4J (núcleo conceptual, genérico sobre `Float` — autodiff-first).
- [ ] HBV-light.
- [ ] Calibración: SCE-UA y DDS.
- [x] Métricas: NSE, KGE (+componentes), logNSE, PBIAS.
- [x] Forzantes: series de precipitación/PET (CSV) vía CLI.
- [ ] Validación split-sample.
- [ ] (v0.2) Semi-distribuido por subcuencas; aporte nival (ver `snowmelt-rs`).

## Arquitectura tentativa
- `rainflow-core`: estados del modelo, integración temporal, optimizadores.
- Targets: native (Rayon para multi-cuenca) + Python (PyO3) + CLI.

## Validación / paridad numérica
Cross-check contra **airGR** (GR4J) y casos CAMELS-CL.

## Venue objetivo
**Environmental Modelling & Software** o **Journal of Hydrology**.

## Conexiones con tu ecosistema
- **Postdoc DICYT**: activo directo para las 15 cuencas BNA.
- Complementa **Hydroflux** (multi-escala: conceptual ↔ físico).
- **Smelt**: emulación/ML de parámetros; **snowmelt-rs**: módulo nival.

## Refinamiento SOTA (2026-06-10)
**Differentiable modeling** es *el* paradigma emergente (δHBV, physics-embedded
learning). Diseñar el core **autodiff-first** para permitir calibración por
gradiente e híbridos física+ML. La investigación del método vive en
`physics-guided-ml`; rainflow es el substrato determinista Rust, acoplado vía
**PyO3**. Alimenta `nowcast` y `snowmelt-rs` como forzantes.

## Próximos pasos al retomar
1. ~~Implementar GR4J + NSE/KGE~~ ✅ (paridad airGR verificada; falta correr una cuenca CAMELS-CL real).
2. Añadir DDS y validar calibración contra airGR (`airGR::Calibration_Michel`).
3. Descargar forzantes CAMELS-CL para 1–2 cuencas BNA y correr split-sample.
4. HBV-light; luego definir formato de subcuencas para el paso semi-distribuido.
