mod forcing;

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Args, Parser};
use rainflow_core::metrics;
use rainflow_core::{Gr4j, Gr4jParams};

#[derive(Parser)]
#[command(
    name = "rainflow",
    version,
    about = "Conceptual rainfall-runoff models in Rust"
)]
enum Cli {
    /// Run GR4J over a CSV forcing file and report goodness-of-fit metrics.
    Run(RunArgs),
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
    }
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
