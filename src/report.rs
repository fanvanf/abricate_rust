//! Report generation — processes BLAST hits and outputs in TSV/CSV/BED/GFF/JSON formats.
//! Implements filtering, deduplication, sorting, and multi-format output.

use crate::blast::{self, BlastHit};
use anyhow::Result;
use serde::Serialize;
use std::io::Write;

/// A processed report row matching the original abricate 14-column format
#[derive(Debug, Clone, Serialize)]
pub struct ReportRow {
    pub file: String,
    pub sequence: String,
    pub start: i64,
    pub end: i64,
    pub strand: String,
    pub gene: String,
    pub coverage: String,
    pub coverage_map: String,
    pub gaps: i64,
    pub perc_coverage: f64,
    pub perc_identity: f64,
    pub database: String,
    pub accession: String,
    pub product: String,
    pub resistance: String,
}

impl ReportRow {
    /// Header labels for TSV/CSV output
    pub fn header() -> Vec<&'static str> {
        vec![
            "FILE",
            "SEQUENCE",
            "START",
            "END",
            "STRAND",
            "GENE",
            "COVERAGE",
            "COVERAGE_MAP",
            "GAPS",
            "%COVERAGE",
            "%IDENTITY",
            "DATABASE",
            "ACCESSION",
            "PRODUCT",
            "RESISTANCE",
        ]
    }

    /// Convert to a vector of strings for tabular output
    pub fn to_fields(&self) -> Vec<String> {
        vec![
            self.file.clone(),
            self.sequence.clone(),
            self.start.to_string(),
            self.end.to_string(),
            self.strand.clone(),
            self.gene.clone(),
            self.coverage.clone(),
            self.coverage_map.clone(),
            self.gaps.to_string(),
            format!("{:.2}", self.perc_coverage),
            format!("{:.2}", self.perc_identity),
            self.database.clone(),
            self.accession.clone(),
            self.product.clone(),
            self.resistance.clone(),
        ]
    }
}

/// Process raw BLAST hits into filtered, deduplicated, sorted report rows
pub fn process_hits(
    hits: Vec<BlastHit>,
    file_name: &str,
    db_name: &str,
    minid: f64,
    mincov: f64,
    nopath: bool,
) -> Vec<ReportRow> {
    let mut rows: Vec<ReportRow> = Vec::new();
    let mut seen: std::collections::HashSet<(String, i64, i64)> = std::collections::HashSet::new();

    // Display file name: optionally strip path
    let display_file = if nopath {
        std::path::Path::new(file_name)
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| file_name.to_string())
    } else {
        file_name.to_string()
    };

    for hit in hits {
        // Calculate coverage
        let pccov = if hit.slen > 0 {
            100.0 * (hit.length - hit.gaps) as f64 / hit.slen as f64
        } else {
            0.0
        };

        // Filter by mincov
        if pccov < mincov {
            continue;
        }

        // Filter by minid (pident)
        if hit.pident < minid {
            continue;
        }

        // Determine strand and start/end
        let (start, end, strand) = if hit.qstart <= hit.qend {
            (hit.qstart, hit.qend, "+")
        } else {
            (hit.qend, hit.qstart, "-")
        };

        // Deduplicate: same qseqid + qstart + qend, keep first
        let key = (hit.qseqid.clone(), start, end);
        if seen.contains(&key) {
            continue;
        }
        seen.insert(key);

        // Parse sseqid (truncated at first space by BLAST) for db/gene/acc/resistance
        // sseqid format: DB~~~GENE~~~ACCESSION~~~RESISTANCE
        let (sdb, gene, acc, resistance) = blast::parse_sseqid(&hit.sseqid);

        // stitle contains the product description (after the ID in BLAST output)
        let product = hit.stitle.clone();

        // Use the actual database name from the hit, or fall back to the provided db_name
        let db_display = if sdb.is_empty() { db_name.to_string() } else { sdb };

        // Coverage string: sstart-send/slen
        let s_start = hit.sstart.min(hit.send);
        let s_end = hit.sstart.max(hit.send);
        let coverage_str = format!("{}-{}/{}", s_start, s_end, hit.slen);

        // Coverage map
        let cov_map = blast::coverage_map(hit.sstart, hit.send, hit.slen);

        rows.push(ReportRow {
            file: display_file.clone(),
            sequence: hit.qseqid,
            start,
            end,
            strand: strand.to_string(),
            gene,
            coverage: coverage_str,
            coverage_map: cov_map,
            gaps: hit.gaps,
            perc_coverage: pccov,
            perc_identity: hit.pident,
            database: db_display,
            accession: acc,
            product,
            resistance,
        });
    }

    // Sort by SEQUENCE name, then START position
    rows.sort_by(|a, b| {
        a.sequence
            .cmp(&b.sequence)
            .then_with(|| a.start.cmp(&b.start))
    });

    rows
}

/// Output format enum
#[derive(Debug, Clone, Copy)]
pub enum OutFormat {
    Tsv,
    Csv,
    Bed,
    Gff,
    Json,
}

impl OutFormat {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s.to_lowercase().as_str() {
            "tsv" => Ok(OutFormat::Tsv),
            "csv" => Ok(OutFormat::Csv),
            "bed" => Ok(OutFormat::Bed),
            "gff" => Ok(OutFormat::Gff),
            "json" => Ok(OutFormat::Json),
            _ => Err(format!("Unknown output format: {}. Use tsv, csv, bed, gff, or json", s)),
        }
    }
}

/// Write report rows in the specified format
pub fn write_report<W: Write>(
    rows: &[ReportRow],
    format: OutFormat,
    noheader: bool,
    writer: &mut W,
) -> Result<()> {
    match format {
        OutFormat::Tsv => write_tsv(rows, noheader, writer)?,
        OutFormat::Csv => write_csv(rows, noheader, writer)?,
        OutFormat::Bed => write_bed(rows, writer)?,
        OutFormat::Gff => write_gff(rows, noheader, writer)?,
        OutFormat::Json => write_json(rows, writer)?,
    }
    Ok(())
}

fn write_tsv<W: Write>(rows: &[ReportRow], noheader: bool, writer: &mut W) -> Result<()> {
    if !noheader {
        let header: Vec<&str> = ReportRow::header();
        writeln!(writer, "#{}", header.join("\t"))?;
    }
    for row in rows {
        let fields = row.to_fields();
        writeln!(writer, "{}", fields.join("\t"))?;
    }
    Ok(())
}

fn write_csv<W: Write>(rows: &[ReportRow], noheader: bool, writer: &mut W) -> Result<()> {
    if !noheader {
        let header: Vec<&str> = ReportRow::header();
        writeln!(writer, "{}", header.join(","))?;
    }
    for row in rows {
        let fields = row.to_fields();
        // Quote fields that contain commas
        let quoted: Vec<String> = fields
            .iter()
            .map(|f| {
                if f.contains(',') || f.contains('"') {
                    format!("\"{}\"", f.replace('"', "\"\""))
                } else {
                    f.clone()
                }
            })
            .collect();
        writeln!(writer, "{}", quoted.join(","))?;
    }
    Ok(())
}

fn write_bed<W: Write>(rows: &[ReportRow], writer: &mut W) -> Result<()> {
    for row in rows {
        // BED format: chrom  start  end  name  score  strand
        // BED uses 0-based start
        let name = if row.gene.is_empty() {
            row.accession.clone()
        } else {
            row.gene.clone()
        };
        let score = format!("{:.1}", row.perc_coverage);
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            row.sequence,
            row.start - 1, // BED is 0-based
            row.end,
            name,
            score,
            row.strand
        )?;
    }
    Ok(())
}

fn write_gff<W: Write>(rows: &[ReportRow], noheader: bool, writer: &mut W) -> Result<()> {
    if !noheader {
        writeln!(writer, "##gff-version 3")?;
    }
    for row in rows {
        // GFF3 format: seqid  source  type  start  end  score  strand  phase  attributes
        let source = "abricate_rust";
        let feature_type = "AMR_gene";
        let score = format!("{:.2}", row.perc_identity);
        let attributes = format!(
            "Name={};Gene={};Database={};Accession={};Coverage={};Identity={};Product={};Resistance={}",
            row.gene,
            row.gene,
            row.database,
            row.accession,
            format!("{:.2}", row.perc_coverage),
            format!("{:.2}", row.perc_identity),
            row.product.replace(';', "\\;"),
            row.resistance.replace(';', "\\;")
        );
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.sequence,
            source,
            feature_type,
            row.start,
            row.end,
            score,
            row.strand,
            ".",
            attributes
        )?;
    }
    Ok(())
}

fn write_json<W: Write>(rows: &[ReportRow], writer: &mut W) -> Result<()> {
    let json = serde_json::to_string_pretty(rows)?;
    writeln!(writer, "{}", json)?;
    Ok(())
}
