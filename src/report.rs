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
    pub gapopen: i64,          // 新增
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
     pub fn to_fields(&self) -> Vec<String> {
        let gaps_str = format!("{}/{}", self.gapopen, self.gaps);
        vec![
            self.file.clone(),
            self.sequence.clone(),
            self.start.to_string(),
            self.end.to_string(),
            self.strand.clone(),
            self.gene.clone(),
            self.coverage.clone(),
            self.coverage_map.clone(),
            gaps_str,                     // 改为 gapopen/gaps
            format!("{:.2}", self.perc_coverage),
            format!("{:.2}", self.perc_identity),
            self.database.clone(),
            self.accession.clone(),
            self.product.clone(),
            self.resistance.clone(),
        ]

    
    }
}
fn format_coverage_map(sstart: i64, send: i64, slen: i64, gapopen: i64) -> String {
    let broken = gapopen > 0;
    let width = if broken { 14 } else { 15 };
    let scale = slen as f64 / width as f64;
    // 计算索引位置，并约束到 [0, width-1]
    let x = ((sstart as f64 / scale).round() as i64)
        .max(0)
        .min((width - 1) as i64) as usize;
    let y = ((send as f64 / scale).round() as i64)
        .max(0)
        .min((width - 1) as i64) as usize;
    let mut map = String::with_capacity(width);
    for i in 0..width {
        let ch = if i >= x && i <= y { '=' } else { '=' };
        map.push(ch);
    }
    if broken {
        let mid = width / 2;
        map.insert(mid, '/');
    }
    map
}

// 修改 process_hits
pub fn process_hits(
    hits: Vec<BlastHit>,
    file_name: &str,
    db_name: &str,
    minid: f64,
    mincov: f64,
    nopath: bool,
) -> Vec<ReportRow> {
    let mut rows = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let display_file = if nopath {
        std::path::Path::new(file_name)
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| file_name.to_string())
    } else {
        file_name.to_string()
    };

    for hit in hits {
        let pccov = if hit.slen > 0 {
            100.0 * (hit.length - hit.gaps) as f64 / hit.slen as f64
        } else {
            0.0
        };
        if pccov < mincov || hit.pident < minid {
            continue;
        }

        // 使用 sstrand 确定方向（假设 hit.sstrand 为 "plus" 或 "minus"）
        let strand = if hit.sstrand == "minus" { "-" } else { "+" };

        // 注意：start/end 仍用 qstart/qend（原始坐标，不交换）
        let (start, end) = (hit.qstart, hit.qend);  // 保留原始值

        let key = (hit.qseqid.clone(), start, end);
        if seen.contains(&key) {
            continue;
        }
        seen.insert(key);

        let (sdb, gene, acc, resistance) = blast::parse_sseqid(&hit.sseqid);

        // 处理 product：若包含 ~~~，去除第一个空格前的标识符
        let product = if hit.stitle.contains("~~~") {
            hit.stitle
                .split_once(' ')
                .map(|(_, desc)| desc.to_string())
                .unwrap_or_else(|| hit.stitle.clone())
        } else {
            hit.stitle.clone()
        };

        let s_start = hit.sstart.min(hit.send);
        let s_end = hit.sstart.max(hit.send);
        let coverage_str = format!("{}-{}/{}", s_start, s_end, hit.slen);
        let cov_map = format_coverage_map(hit.sstart, hit.send, hit.slen, hit.gapopen);

        let db_display = if sdb.is_empty() { db_name.to_string() } else { sdb };

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
            gapopen: hit.gapopen,        // 新增
            perc_coverage: pccov,
            perc_identity: hit.pident,
            database: db_display,
            accession: acc,
            product,
            resistance,
        });
    }

    rows.sort_by(|a, b| a.sequence.cmp(&b.sequence).then_with(|| a.start.cmp(&b.start)));
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
