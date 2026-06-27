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
- [x] HBV-light (rutina nival grado-día con temperatura opcional; en las cuencas
      CAMELS-CL pluviales supera a GR4J: NSE val 0.73–0.77 vs 0.64–0.74).
- [x] Calibración: DDS (validada vs `airGR::Calibration_Michel`: NSE 0.7956 vs 0.7957)
      y SCE-UA (Duan et al. 1992); seleccionables con `--algorithm dds|sce`.
      **Desde 2026-06: DDS/SCE-UA se consumen del substrato `forge` (`forge-core`,
      path dep); `calibrate.rs` solo conserva la API y delega los optimizadores.
      API pública intacta, paridad airGR y RNG bit-idénticos.**
- [x] Métricas: NSE, KGE (+componentes), logNSE, PBIAS.
- [x] Forzantes: series de precipitación/PET (CSV) vía CLI.
- [x] Validación split-sample (`rainflow split-sample`) sobre 2 cuencas CAMELS-CL
      pluviales casi naturales (8123001 Itata en Cholguán, 7330001 Perquilauquén
      en San Manuel): KGE validación 0.76–0.82. Ver `data/camels-cl/README.md`.
- [ ] (v0.2) Semi-distribuido por subcuencas; aporte nival (ver `snowmelt-rs`).

## Arquitectura tentativa
- `rainflow-core`: estados del modelo, integración temporal; optimizadores
  (DDS/SCE-UA) vía `forge-core`.
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
1. ~~Implementar GR4J + NSE/KGE~~ ✅ (paridad airGR verificada).
2. ~~Añadir DDS y validar calibración contra airGR~~ ✅ (mismo óptimo que Calibration_Michel).
3. ~~CAMELS-CL + split-sample~~ ✅ (2 cuencas BNA pluviales, KGE val 0.76–0.82).
4. ~~HBV-light~~ ✅ (sin paridad externa: HBV-light es software GUI; validado por
   invariantes —balance de masa exacto, nieve, cotas— y benchmark vs GR4J).
5. ~~Cuenca nival CAMELS-CL~~ ✅ (4511002 Las Ramadas + 4703002 Cuncumén: rutina
   nival sube NSE val de ≤0.23 a 0.31–0.62; TT/SFCF calibran alto por temperatura
   agregada en cuencas de alto relieve → motiva bandas de elevación en v0.2).
6. ~~Bandas de elevación / semi-distribuido~~ ✅ (`--model hbv-bands`: nieve+suelo
   por banda con TCALT/PCALT, respuesta y ruteo compartidos; 1 banda = modelo
   agregado exacto. En Las Ramadas TT vuelve a ~0°C físico y NSE val 0.33→0.51).
7. ~~SCE-UA~~ ✅ (Duan et al. 1992; concuerda con DDS en GR4J salvo redondeo).
8. ~~Calibrar TCALT/PCALT en bandas~~ ✅ (arregló 4703002: val NSE 0.23→0.76 con
   lapse calibrado; bandas ahora superan al agregado de forma robusta).
9. ~~Geometría hipsométrica de bandas~~ ✅ (`--hypsometry "min,median,max"`:
   bandas equi-área con elevación leída de la curva hipsométrica; iguala/supera
   a las bandas a-ojo, +0.05–0.11 NSE val en 4511002, sin arbitrariedad).
   El core `equal_area_from_hypsometry` acepta también curva densa de DEM.
10. ~~Curva hipsométrica desde DEM~~ ✅ (Copernicus GLO-30 recortado por polígono
    CAMELS-CL oficial vía gateway GIS; `--hypsometry-file`; el recorte reproduce
    los atributos CAMELS-CL a pocos metros. Skill comparable a la reconstrucción
    de cuantiles —el lapse calibrado absorbe el error de geometría— pero objetiva).
11. ~~CI + LICENSE + PyO3~~ ✅.
12. ~~Routing por subcuencas~~ ✅ (módulo `routing`: Muskingum + RiverNetwork
    árbol de drenaje; validado en Itata anidado Cholguán→Balsa Nueva Aldea:
    distribuido 2-subcuencas supera al lumped +0.06–0.08 NSE val. K~0.1d: a
    paso diario el tránsito es casi pasante, la ganancia es la desagregación
    espacial de forzante/parámetros, no el Muskingum per se).
13. ~~Autodiff (ADN del proyecto)~~ ✅ (crate `autodiff` forward-mode implementa
    num_traits::Float → pasa por GR4J sin tocar el core; gradiente analítico
    valida vs diferencias finitas ~1e-9; Adam recupera params verdaderos exacto.
    `examples/gradient_calibration.rs`. Substrato para híbridos física+ML δHBV).
14. ~~Warm-start desde estado~~ ✅ (`run_from` en GR4J/HBV/HbvBands: corre desde
    un estado dado, continuidad exacta 1e-12; PyO3 expone Gr4jState/HbvState
    opacos con semántica de clonado → nowcast asienta estado, snapshot, fan-out
    de escenarios). rainflow listo como motor de nowcast.
15. ~~Benchmarks (criterion)~~ ✅ (benches/models.rs + calibration.rs. GR4J ~1.5ms/27yr,
    HBV ~0.68ms (GR4J paga tanh+2 powf/paso), routing ~78µs. Calibración ~5s/cuenca
    → 15 cuencas BNA serial ~75s, motiva Rayon. Tabla en README. Hallazgo: cada eval
    de calibración aloca qsim+pairs; reuso de buffers es optimización futura).
16. ~~Rayon multi-cuenca~~ ✅ (`rainflow batch`: calibra N cuencas en paralelo;
    core puro → rayon vive solo en el CLI. Speedup 3.86× en 4 cuencas / 16 threads,
    serial 2m51s → 44s. Path operacional para las 15 cuencas BNA).
17. ~~serde-state~~ ✅ (feature `serde` opcional en rainflow-core: deriva Serialize/
    Deserialize en params/estados/config; PyO3 expone to_json/from_json en
    Gr4jState/HbvState. Round-trip exacto desde disco verificado en Rust y Python.
    CI corre `cargo test --features serde`. Core liviano por defecto).
18. Próximo: publicar crates.io/PyPI; paper EMS/JH; arrancar nowcast vía PyO3.
5. Caso interesante para el paper: 8123001 muestra equifinalidad + no-estacionariedad
   (megasequía post-2010) — benchmark para calibración por gradiente/regularizada.
