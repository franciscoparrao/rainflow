//! CSV forcing reader. Expected columns (case-insensitive, flexible names):
//! date, precipitation, PET and optionally observed discharge.

use std::path::Path;

use anyhow::{Context, Result, bail};

#[derive(Debug)]
pub struct Forcing {
    pub dates: Vec<String>,
    pub precip: Vec<f64>,
    pub pet: Vec<f64>,
    pub qobs: Option<Vec<f64>>,
}

const DATE_NAMES: &[&str] = &["date", "fecha", "time", "dates"];
const PRECIP_NAMES: &[&str] = &[
    "p",
    "precip",
    "precipitation",
    "pr",
    "pcp",
    "prcp",
    "precip_mm",
];
const PET_NAMES: &[&str] = &["pet", "etp", "e", "evap", "pet_mm", "pe"];
const QOBS_NAMES: &[&str] = &[
    "q",
    "qobs",
    "q_mm",
    "qmm",
    "discharge",
    "caudal",
    "streamflow",
];

fn find_column(headers: &[String], names: &[&str]) -> Option<usize> {
    headers
        .iter()
        .position(|h| names.contains(&h.trim().to_lowercase().as_str()))
}

fn parse_value(raw: &str) -> Result<f64> {
    let v = raw.trim();
    if v.is_empty() || v.eq_ignore_ascii_case("na") || v.eq_ignore_ascii_case("nan") {
        return Ok(f64::NAN);
    }
    v.parse::<f64>()
        .with_context(|| format!("invalid number: {v:?}"))
}

pub fn read_csv(path: &Path) -> Result<Forcing> {
    let mut reader = csv::Reader::from_path(path)
        .with_context(|| format!("cannot open forcing file {}", path.display()))?;
    let headers: Vec<String> = reader
        .headers()
        .context("forcing file has no header row")?
        .iter()
        .map(str::to_owned)
        .collect();

    let date_col = find_column(&headers, DATE_NAMES);
    let Some(p_col) = find_column(&headers, PRECIP_NAMES) else {
        bail!("no precipitation column found (tried {PRECIP_NAMES:?}) in {headers:?}");
    };
    let Some(pet_col) = find_column(&headers, PET_NAMES) else {
        bail!("no PET column found (tried {PET_NAMES:?}) in {headers:?}");
    };
    let q_col = find_column(&headers, QOBS_NAMES);

    let mut forcing = Forcing {
        dates: Vec::new(),
        precip: Vec::new(),
        pet: Vec::new(),
        qobs: q_col.map(|_| Vec::new()),
    };

    for (i, record) in reader.records().enumerate() {
        let record = record.with_context(|| format!("error reading row {}", i + 2))?;
        let row_ctx = || format!("row {} of {}", i + 2, path.display());

        forcing.dates.push(
            date_col
                .and_then(|c| record.get(c))
                .map(str::to_owned)
                .unwrap_or_else(|| (i + 1).to_string()),
        );
        let p = parse_value(record.get(p_col).unwrap_or("")).with_context(row_ctx)?;
        let pet = parse_value(record.get(pet_col).unwrap_or("")).with_context(row_ctx)?;
        if !p.is_finite() || !pet.is_finite() {
            bail!(
                "missing precipitation/PET at {} — gaps are only allowed in qobs",
                row_ctx()
            );
        }
        forcing.precip.push(p);
        forcing.pet.push(pet);
        if let (Some(c), Some(qs)) = (q_col, forcing.qobs.as_mut()) {
            qs.push(parse_value(record.get(c).unwrap_or("")).with_context(row_ctx)?);
        }
    }

    if forcing.precip.is_empty() {
        bail!("forcing file {} has no data rows", path.display());
    }
    Ok(forcing)
}
