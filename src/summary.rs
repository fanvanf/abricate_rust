//! Summary mode — reads multiple abricate TSV reports and generates a gene
//! presence/absence matrix (rows = files, columns = genes).

use anyhow::{Context, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::path::Path;

/// Parse an abricate TSV report file and extract (gene, value) pairs
fn parse_abricate_report(path: &str, use_identity: bool) -> Result<(String, Vec<(String, f64)>)> {
    let content = if path == "-" {
        use std::io::Read;
        let mut s = String::new();
        std::io::stdin().read_to_string(&mut s)?;
        s
    } else {
        std::fs::read_to_string(path)
            .with_context(|| format!("Cannot read report file: {}", path))?
    };

    // Extract display name (strip path and extension)
    let display_name = if path == "-" {
        "stdin".to_string()
    } else {
        Path::new(path)
            .file_stem()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| path.to_string())
    };

    let mut genes = Vec::new();

    for line in content.lines() {
        // Skip header lines (start with #) and empty lines
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 15 {
            continue;
        }

        // Column indices (0-based):
        // 0=FILE, 1=SEQUENCE, 2=START, 3=END, 4=STRAND, 5=GENE,
        // 6=COVERAGE, 7=COVERAGE_MAP, 8=GAPS, 9=%COVERAGE, 10=%IDENTITY,
        // 11=DATABASE, 12=ACCESSION, 13=PRODUCT, 14=RESISTANCE
        let gene = fields[5].to_string();
        let value = if use_identity {
            fields[10].parse::<f64>().unwrap_or(0.0)
        } else {
            fields[9].parse::<f64>().unwrap_or(0.0)
        };

        if !gene.is_empty() {
            genes.push((gene, value));
        }
    }

    Ok((display_name, genes))
}

/// Generate summary matrix from multiple abricate report files
pub fn generate_summary<W: Write>(
    files: &[String],
    use_identity: bool,
    noheader: bool,
    writer: &mut W,
) -> Result<()> {
    // Collect all data: file_name -> {gene -> value}
    let mut all_data: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();
    let mut all_genes: BTreeSet<String> = BTreeSet::new();

    for file in files {
        let (name, genes) = parse_abricate_report(file, use_identity)?;

        let mut gene_map: BTreeMap<String, f64> = BTreeMap::new();
        for (gene, value) in genes {
            // Keep the maximum value if gene appears multiple times
            let entry = gene_map.entry(gene.clone()).or_insert(0.0);
            if value > *entry {
                *entry = value;
            }
            all_genes.insert(gene);
        }

        all_data.insert(name, gene_map);
    }

    // Sort file names
    let sorted_files: Vec<&String> = all_data.keys().collect();
    let sorted_genes: Vec<&String> = all_genes.iter().collect();

    // Write header
    if !noheader {
        let mut header = vec!["#FILE".to_string()];
        for g in &sorted_genes {
            header.push((**g).clone());
        }
        writeln!(writer, "{}", header.join("\t"))?;
    }

    // Write data rows
    for file_name in &sorted_files {
        let gene_map = &all_data[*file_name];
        let mut row = vec![(*file_name).clone()];
        for g in &sorted_genes {
            if let Some(val) = gene_map.get(*g) {
                row.push(format!("{:.1}", val));
            } else {
                row.push(".".to_string());
            }
        }
        writeln!(writer, "{}", row.join("\t"))?;
    }

    Ok(())
}
