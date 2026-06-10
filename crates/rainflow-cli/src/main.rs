mod forcing;

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Args, Parser};
use rainflow_core::calibrate::{self, DdsConfig, Objective};
use rainflow_core::metrics;
use rainflow_core::{Gr4j, Gr4jParams, Hbv};

#[derive(Parser)]
#[command(
    name = "rainflow",
    version,
    about = "Conceptual rainfall-runoff models in Rust"
)]
enum Cli {
    /// Run GR4J over a CSV forcing file and report goodness-of-fit metrics.
    Run(RunArgs),
    /// Calibrate GR4J with DDS against observed discharge.
    Calibrate(CalibrateArgs),
    /// Split-sample test: calibrate on each half of the record, validate on
    /// the other (Klemeš 1986).
    SplitSample(CalibrateArgs),
}

#[derive(Args)]
struct CalibrateArgs {
    /// CSV with columns: date, precipitation, PET, qobs (mm) and optionally
    /// temperature (°C, enables the HBV snow routine)
    #[arg(long)]
    forcing: PathBuf,

    /// Model to calibrate
    #[arg(long, value_parser = ["gr4j", "hbv"], default_value = "gr4j")]
    model: String,

    /// Objective to maximize
    #[arg(long, value_parser = ["nse", "kge", "lognse"], default_value = "nse")]
    objective: String,

    /// DDS evaluation budget
    #[arg(long, default_value_t = 2000)]
    iterations: usize,

    /// RNG seed (same seed => same result)
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Warm-up steps excluded from the objective
    #[arg(long, default_value_t = 365)]
    warmup: usize,
}

#[derive(Args)]
struct RunArgs {
    /// CSV with columns: date, precipitation, PET [, qobs] (mm)
    #[arg(long)]
    forcing: PathBuf,

    /// GR4J x1 — production store capacity (mm)
    #[arg(long)]
    x1: f64,
    /// GR4J x2 — groundwater exchange coefficient (mm/day)
    #[arg(long, allow_hyphen_values = true)]
    x2: f64,
    /// GR4J x3 — routing store capacity (mm)
    #[arg(long)]
    x3: f64,
    /// GR4J x4 — unit hydrograph time base (days)
    #[arg(long)]
    x4: f64,

    /// Warm-up steps excluded from the metrics
    #[arg(long, default_value_t = 365)]
    warmup: usize,

    /// Write simulated discharge (date,qsim) to this CSV
    #[arg(long)]
    output: Option<PathBuf>,
}

fn main() -> Result<()> {
    match Cli::parse() {
        Cli::Run(args) => run(args),
        Cli::Calibrate(args) => calibrate(args),
        Cli::SplitSample(args) => split_sample(args),
    }
}

fn split_sample(args: CalibrateArgs) -> Result<()> {
    let forcing = forcing::read_csv(&args.forcing)?;
    let qobs = forcing
        .qobs
        .as_deref()
        .context("split-sample requires an observed discharge column (qobs)")?;
    let n = qobs.len();
    anyhow::ensure!(
        args.warmup < n / 2,
        "warm-up ({}) must be shorter than half the series ({})",
        args.warmup,
        n / 2
    );

    let objective = parse_objective(&args.objective)?;
    let config = DdsConfig {
        max_iter: args.iterations,
        seed: args.seed,
        ..Default::default()
    };
    let mid = n / 2;

    // Mask observations outside a period with NaN; the metrics skip them.
    let mask = |range: std::ops::Range<usize>| -> Vec<f64> {
        qobs.iter()
            .enumerate()
            .map(|(i, &v)| if range.contains(&i) { v } else { f64::NAN })
            .collect()
    };
    let halves = [
        (
            "A (first half)",
            "B (second half)",
            mask(args.warmup..mid),
            mask(mid..n),
        ),
        (
            "B (second half)",
            "A (first half)",
            mask(mid..n),
            mask(args.warmup..mid),
        ),
    ];

    println!(
        "Split-sample test of {} ({} steps, split at {mid}, objective {}, {} evaluations)",
        args.model, n, args.objective, args.iterations
    );
    for (cal_name, val_name, cal_obs, val_obs) in &halves {
        let cal = calibrate_dispatch(
            &args.model,
            &forcing,
            cal_obs,
            args.warmup,
            objective,
            &config,
        )?;
        let val = objective_value(objective, &val_obs[args.warmup..], &cal.qsim[args.warmup..]);
        println!(
            "  cal {cal_name}: {:.4}  ->  val {val_name}: {}   [{}]",
            cal.value,
            val.map_or_else(|| "n/a".into(), |v| format!("{v:.4}")),
            cal.desc
        );
    }
    Ok(())
}

/// Result of calibrating either model: optimum, its description and the
/// simulated discharge at the optimum over the full series.
struct CalOutcome {
    value: f64,
    evaluations: usize,
    desc: String,
    qsim: Vec<f64>,
}

fn calibrate_dispatch(
    model: &str,
    forcing: &forcing::Forcing,
    qobs: &[f64],
    warmup: usize,
    objective: Objective,
    config: &DdsConfig,
) -> Result<CalOutcome> {
    match model {
        "gr4j" => {
            let cal = calibrate::calibrate_gr4j(
                &forcing.precip,
                &forcing.pet,
                qobs,
                warmup,
                objective,
                &calibrate::gr4j_default_bounds(),
                config,
            )?;
            let qsim = Gr4j::new(cal.params)?.run(&forcing.precip, &forcing.pet)?;
            let p = cal.params;
            Ok(CalOutcome {
                value: cal.value,
                evaluations: cal.evaluations,
                desc: format!(
                    "x1={:.1} x2={:.2} x3={:.1} x4={:.2}",
                    p.x1, p.x2, p.x3, p.x4
                ),
                qsim,
            })
        }
        "hbv" => {
            let temp = forcing.temp.as_deref();
            let cal = calibrate::calibrate_hbv(
                &forcing.precip,
                &forcing.pet,
                temp,
                qobs,
                warmup,
                objective,
                config,
            )?;
            let qsim = Hbv::new(cal.params)?.run(&forcing.precip, &forcing.pet, temp)?;
            let p = cal.params;
            let snow = if temp.is_some() {
                format!("tt={:.2} cfmax={:.2} sfcf={:.2} ", p.tt, p.cfmax, p.sfcf)
            } else {
                String::new()
            };
            Ok(CalOutcome {
                value: cal.value,
                evaluations: cal.evaluations,
                desc: format!(
                    "{snow}fc={:.0} lp={:.2} beta={:.2} k0={:.3} k1={:.3} k2={:.4} \
                     uzl={:.1} perc={:.2} maxbas={:.2}",
                    p.fc, p.lp, p.beta, p.k0, p.k1, p.k2, p.uzl, p.perc, p.maxbas
                ),
                qsim,
            })
        }
        other => anyhow::bail!("unknown model {other:?}"),
    }
}

fn parse_objective(name: &str) -> Result<Objective> {
    match name {
        "nse" => Ok(Objective::Nse),
        "kge" => Ok(Objective::Kge),
        "lognse" => Ok(Objective::LogNse),
        other => anyhow::bail!("unknown objective {other:?}"),
    }
}

fn objective_value(objective: Objective, obs: &[f64], sim: &[f64]) -> Option<f64> {
    match objective {
        Objective::Nse => metrics::nse(obs, sim).ok(),
        Objective::Kge => metrics::kge(obs, sim).ok(),
        Objective::LogNse => metrics::log_nse(obs, sim).ok(),
    }
}

fn calibrate(args: CalibrateArgs) -> Result<()> {
    let forcing = forcing::read_csv(&args.forcing)?;
    let qobs = forcing
        .qobs
        .as_deref()
        .context("calibration requires an observed discharge column (qobs)")?;

    let objective = parse_objective(&args.objective)?;
    let config = DdsConfig {
        max_iter: args.iterations,
        seed: args.seed,
        ..Default::default()
    };

    let cal = calibrate_dispatch(&args.model, &forcing, qobs, args.warmup, objective, &config)?;

    println!(
        "DDS calibration of {} ({} evaluations, seed {}, warm-up {})",
        args.model, cal.evaluations, args.seed, args.warmup
    );
    println!("  best params: {}", cal.desc);
    println!("  best {:>6}: {:.4}", args.objective, cal.value);

    // Report the full metric suite at the optimum.
    let obs = &qobs[args.warmup..];
    let sim = &cal.qsim[args.warmup..];
    print_metric("NSE", metrics::nse(obs, sim));
    print_metric("KGE", metrics::kge(obs, sim));
    print_metric("logNSE", metrics::log_nse(obs, sim));
    print_metric("PBIAS%", metrics::pbias(obs, sim));
    Ok(())
}

fn run(args: RunArgs) -> Result<()> {
    let forcing = forcing::read_csv(&args.forcing)?;
    let n = forcing.precip.len();
    anyhow::ensure!(
        args.warmup < n,
        "warm-up ({}) must be shorter than the series ({n} steps)",
        args.warmup
    );

    let model = Gr4j::new(Gr4jParams {
        x1: args.x1,
        x2: args.x2,
        x3: args.x3,
        x4: args.x4,
    })?;
    let qsim = model.run(&forcing.precip, &forcing.pet)?;

    println!(
        "GR4J  x1={} x2={} x3={} x4={}  |  {n} steps, warm-up {}",
        args.x1, args.x2, args.x3, args.x4, args.warmup
    );

    if let Some(qobs) = &forcing.qobs {
        let obs = &qobs[args.warmup..];
        let sim = &qsim[args.warmup..];
        print_metric("NSE", metrics::nse(obs, sim));
        print_metric("KGE", metrics::kge(obs, sim));
        print_metric("logNSE", metrics::log_nse(obs, sim));
        print_metric("PBIAS%", metrics::pbias(obs, sim));
        if let Ok(c) = metrics::kge_components(obs, sim) {
            println!(
                "  KGE components: r={:.4} alpha={:.4} beta={:.4}",
                c.r, c.alpha, c.beta
            );
        }
    } else {
        println!("  (no observed discharge column — metrics skipped)");
    }

    if let Some(path) = &args.output {
        let mut writer = csv::Writer::from_path(path)
            .with_context(|| format!("cannot write {}", path.display()))?;
        writer.write_record(["date", "qsim_mm"])?;
        for (date, q) in forcing.dates.iter().zip(&qsim) {
            writer.write_record([date.as_str(), &format!("{q:.6}")])?;
        }
        writer.flush()?;
        println!("  simulated discharge written to {}", path.display());
    }

    Ok(())
}

fn print_metric(name: &str, value: Result<f64, metrics::MetricsError>) {
    match value {
        Ok(v) => println!("  {name:<7} {v:>8.4}"),
        Err(e) => println!("  {name:<7} n/a ({e})"),
    }
}
