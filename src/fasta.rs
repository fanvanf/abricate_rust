//! Integrated any2fasta logic — converts various sequence formats to uppercase FASTA.
//! Refactored from the standalone any2fasta tool to write to an internal buffer
//! instead of stdout, enabling in-memory pipeline without external process calls.

use anyhow::{Context, Result};
use std::cell::Cell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};

// ---------------------------------------------------------------------------
// Open input with automatic decompression
// ---------------------------------------------------------------------------
fn open_input(path: &str) -> Result<Box<dyn Read>> {
    if path == "-" {
        return Ok(Box::new(io::stdin()));
    }
    let file = File::open(path).with_context(|| format!("Cannot open '{}'", path))?;
    let path_lower = path.to_lowercase();

    if path_lower.ends_with(".gz") {
        Ok(Box::new(flate2::read::MultiGzDecoder::new(file)))
    } else if path_lower.ends_with(".bz2") {
        Ok(Box::new(bzip2::read::BzDecoder::new(file)))
    } else if path_lower.ends_with(".zst") || path_lower.ends_with(".zstd") {
        Ok(Box::new(zstd::stream::read::Decoder::new(file)?))
    } else if path_lower.ends_with(".zip") {
        let mut archive = zip::ZipArchive::new(file)?;
        let mut entry = archive.by_index(0).context("Zip archive has no entries")?;
        let mut content = Vec::new();
        entry.read_to_end(&mut content)?;
        Ok(Box::new(io::Cursor::new(content)))
    } else {
        Ok(Box::new(file))
    }
}

// ---------------------------------------------------------------------------
// Internal config for the converter
// ---------------------------------------------------------------------------
#[allow(dead_code)]
struct ConvConfig {
    quiet: bool,
    uppercase: bool,
    protein: Cell<bool>,
}

impl ConvConfig {
    fn new() -> Self {
        ConvConfig {
            quiet: true,
            uppercase: true,
            protein: Cell::new(false),
        }
    }
}

fn is_protein(seq: &str) -> bool {
    let len = seq.len();
    if len == 0 {
        return false;
    }
    let non_atgc = seq
        .chars()
        .filter(|c| !"ATCGatcg".contains(*c))
        .count() as f64;
    non_atgc / (len as f64) > 0.666
}

fn purify_seq(config: &ConvConfig, seq: &str) -> String {
    let mut s = seq.to_string();
    if config.uppercase {
        s = s.to_uppercase();
    }
    s
}

fn purify_id(seq: &str) -> String {
    // For abricate use, strip descriptions (take first word of header)
    if let Some(p) = seq.find(|c: char| c == ' ' || c == '\t') {
        seq[..p].to_string()
    } else {
        seq.to_string()
    }
}

// ---------------------------------------------------------------------------
// Parser functions — write to a Writer
// ---------------------------------------------------------------------------
fn parse_fasta<W: Write>(config: &ConvConfig, lines: &mut Vec<String>, out: &mut W) -> Result<usize> {
    let mut count = 0;
    let mut seq_line_num = 0;
    let mut first_seq = true;
    for line in lines.iter_mut() {
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('>') {
            *line = purify_id(line);
            count += 1;
            seq_line_num = 0;
            writeln!(out, "{}", line)?;
        } else {
            seq_line_num += 1;
            if count == 1 && seq_line_num == 1 && first_seq {
                let is_prot = is_protein(line);
                config.protein.set(is_prot);
                first_seq = false;
            }
            let pur = purify_seq(config, line);
            writeln!(out, "{}", pur)?;
        }
    }
    Ok(count)
}

fn parse_fastq<W: Write>(config: &ConvConfig, lines: &mut Vec<String>, out: &mut W) -> Result<usize> {
    let mut count = 0;
    let mut i = 0;
    while i + 3 < lines.len() {
        let seq_line = lines[i + 1].clone();
        let id_line = purify_id(&lines[i]);
        writeln!(out, ">{}", &id_line[1..])?;
        let pur = purify_seq(config, &seq_line);
        writeln!(out, "{}", pur)?;
        count += 1;
        i += 4;
    }
    Ok(count)
}

fn parse_gff<W: Write>(config: &ConvConfig, lines: &mut Vec<String>, out: &mut W) -> Result<usize> {
    let mut fasta_lines = Vec::new();
    let mut in_fasta = false;
    for line in lines.iter() {
        if line.starts_with("##FASTA") {
            in_fasta = true;
            continue;
        }
        if in_fasta && (line.starts_with('>') || !line.trim().is_empty()) {
            fasta_lines.push(line.clone());
        }
    }
    if fasta_lines.is_empty() {
        return Ok(0);
    }
    parse_fasta(config, &mut fasta_lines, out)
}

fn parse_gfa<W: Write>(config: &ConvConfig, lines: &mut Vec<String>, out: &mut W) -> Result<usize> {
    let mut count = 0;
    for line in lines {
        if line.starts_with('S') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                writeln!(out, ">{}", parts[1])?;
                let pur = purify_seq(config, parts[2]);
                writeln!(out, "{}", pur)?;
                count += 1;
            }
        }
    }
    Ok(count)
}

fn parse_genbank<W: Write>(config: &ConvConfig, lines: &mut Vec<String>, out: &mut W) -> Result<usize> {
    let mut acc = String::new();
    let mut dna = String::new();
    let mut in_seq = false;
    let mut count = 0;

    for line in lines {
        let line = line.trim_end();
        if line.starts_with("//") {
            if !acc.is_empty() && !dna.is_empty() {
                writeln!(out, ">{}", acc)?;
                let pur = purify_seq(config, &dna);
                writeln!(out, "{}", pur)?;
                count += 1;
            }
            in_seq = false;
            dna.clear();
            acc.clear();
            continue;
        }
        if line.starts_with("ORIGIN") {
            in_seq = true;
            continue;
        }
        if in_seq {
            if line.len() > 10 {
                let s = &line[10..];
                let s = s.replace(' ', "").replace('\t', "");
                dna.push_str(&s);
            }
        } else if line.starts_with("LOCUS") {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() >= 2 {
                acc = fields[1].to_string();
                if line.contains(" aa ") || line.contains(" AA ") {
                    config.protein.set(true);
                }
            }
        }
    }
    if !acc.is_empty() && !dna.is_empty() {
        writeln!(out, ">{}", acc)?;
        let pur = purify_seq(config, &dna);
        write!(out, "{}", pur)?;
        count += 1;
    }
    Ok(count)
}

fn parse_embl<W: Write>(config: &ConvConfig, lines: &mut Vec<String>, out: &mut W) -> Result<usize> {
    let mut acc = String::new();
    let mut dna = String::new();
    let mut in_seq = false;
    let mut count = 0;

    for line in lines {
        let line = line.trim_end();
        if line.starts_with("//") {
            if !acc.is_empty() && !dna.is_empty() {
                writeln!(out, ">{}", acc)?;
                let pur = purify_seq(config, &dna);
                writeln!(out, "{}", pur)?;
                count += 1;
            }
            in_seq = false;
            dna.clear();
            acc.clear();
            continue;
        }
        if line.starts_with("SQ") {
            in_seq = true;
            if line.contains(" AA ") || line.contains(" aa ") {
                config.protein.set(true);
            }
            continue;
        }
        if in_seq {
            let mut s = line.to_string();
            s.retain(|c| !c.is_whitespace() && !c.is_ascii_digit());
            if !s.is_empty() {
                dna.push_str(&s);
            }
        } else if line.starts_with("ID") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                let id_part = parts[1];
                let id = id_part.trim_end_matches(';');
                acc = id.to_string();
            }
        }
    }
    if !acc.is_empty() && !dna.is_empty() {
        writeln!(out, ">{}", acc)?;
        let pur = purify_seq(config, &dna);
        write!(out, "{}", pur)?;
        count += 1;
    }
    Ok(count)
}

fn parse_clustal<W: Write>(config: &ConvConfig, lines: &mut Vec<String>, out: &mut W) -> Result<usize> {
    let mut seqs: HashMap<String, String> = HashMap::new();
    let mut order = Vec::new();

    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with("CLUSTAL") || trimmed.starts_with("MUSCLE") {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 2 {
            continue;
        }
        let id = fields[0];
        let seq = fields[1];
        let seq_clean = seq.trim_end_matches(|c: char| c.is_ascii_digit());
        if seq_clean.is_empty() {
            continue;
        }
        if !seq_clean.chars().all(|c| c.is_ascii_alphabetic() || c == '-') {
            continue;
        }
        if !seqs.contains_key(id) {
            order.push(id.to_string());
        }
        let entry = seqs.entry(id.to_string()).or_default();
        entry.push_str(seq_clean);
    }

    for id in &order {
        if let Some(seq) = seqs.get(id) {
            writeln!(out, ">{}", id)?;
            let pur = purify_seq(config, seq);
            writeln!(out, "{}", pur)?;
        }
    }
    Ok(order.len())
}

fn parse_stockholm<W: Write>(config: &ConvConfig, lines: &mut Vec<String>, out: &mut W) -> Result<usize> {
    let mut seqs: HashMap<String, String> = HashMap::new();
    let mut order = Vec::new();

    for line in lines {
        let trimmed = line.trim();
        if trimmed.starts_with('#') || trimmed.starts_with("//") || trimmed.is_empty() {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 2 {
            continue;
        }
        let id = fields[0];
        let seq = fields[1];
        let seq_clean = seq.trim_end_matches(|c: char| c.is_ascii_digit());
        if seq_clean.is_empty() {
            continue;
        }
        let seq_replaced = seq_clean.replace('.', "-");
        if !seq_replaced.chars().all(|c| c.is_ascii_alphabetic() || c == '-') {
            continue;
        }
        if !seqs.contains_key(id) {
            order.push(id.to_string());
        }
        let entry = seqs.entry(id.to_string()).or_default();
        entry.push_str(&seq_replaced);
    }

    for id in &order {
        if let Some(seq) = seqs.get(id) {
            writeln!(out, ">{}", id)?;
            let pur = purify_seq(config, seq);
            write!(out, "{}", pur)?;
        }
    }
    Ok(order.len())
}

fn parse_pdb<W: Write>(config: &ConvConfig, lines: &mut Vec<String>, out: &mut W) -> Result<usize> {
    let aa_map: HashMap<&str, char> = [
        ("ALA", 'A'), ("ARG", 'R'), ("ASN", 'N'), ("ASP", 'D'),
        ("CYS", 'C'), ("GLU", 'E'), ("GLN", 'Q'), ("GLY", 'G'),
        ("HIS", 'H'), ("ILE", 'I'), ("LEU", 'L'), ("LYS", 'K'),
        ("MET", 'M'), ("PHE", 'F'), ("PRO", 'P'), ("SER", 'S'),
        ("THR", 'T'), ("TRP", 'W'), ("TYR", 'Y'), ("VAL", 'V'),
        ("SEC", 'U'), ("PYL", 'O'), ("ASX", 'B'), ("GLX", 'Z'),
        ("XAA", 'X'), ("TER", '*'),
    ]
    .iter()
    .cloned()
    .collect();

    let mut seqs: HashMap<String, String> = HashMap::new();
    let mut pdb_id = "unknown".to_string();

    for line in &*lines {
        if line.starts_with("HEADER") {
            if let Some(id) = line.split_whitespace().last() {
                pdb_id = id.to_string();
            }
            break;
        }
    }

    for line in &*lines {
        if line.starts_with("SEQRES") {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 4 {
                continue;
            }
            let chain = fields[2];
            let key = format!("{}-{}", pdb_id, chain);
            let seq = seqs.entry(key).or_default();
            for &aa3 in &fields[4..] {
                seq.push(*aa_map.get(aa3).unwrap_or(&'X'));
            }
        }
    }

    for (id, seq) in &seqs {
        writeln!(out, ">{}", id)?;
        let pur = purify_seq(config, seq);
        writeln!(out, "{}", pur)?;
    }
    Ok(seqs.len())
}

// ---------------------------------------------------------------------------
// Format detection
// ---------------------------------------------------------------------------
type ParserFn<W> = fn(&ConvConfig, &mut Vec<String>, &mut W) -> Result<usize>;

#[allow(dead_code)]
struct Format<W: Write> {
    name: &'static str,
    detect: fn(&str) -> bool,
    parser: ParserFn<W>,
}

fn detect_genbank(line: &str) -> bool {
    line.starts_with("LOCUS ")
}
fn detect_uniprot(line: &str) -> bool {
    line.starts_with("ID") && line.contains("Reviewed")
}
fn detect_embl(line: &str) -> bool {
    line.starts_with("ID ") && !line.contains("Reviewed")
}
fn detect_fasta(line: &str) -> bool {
    line.starts_with('>')
}
fn detect_fastq(line: &str) -> bool {
    line.starts_with('@')
}
fn detect_gff(line: &str) -> bool {
    line.starts_with("##gff")
}
fn detect_clustal(line: &str) -> bool {
    line.starts_with("CLUST") || line.starts_with("MUSCL")
}
fn detect_stockholm(line: &str) -> bool {
    line.starts_with("# STOCKHOLM ")
}
fn detect_gfa(line: &str) -> bool {
    if line.len() < 2 {
        return false;
    }
    let first = line.chars().next().unwrap();
    let second = line.chars().nth(1).unwrap();
    first.is_ascii_uppercase() && (second == '\t' || second == ' ')
}
fn detect_pdb(line: &str) -> bool {
    line.starts_with("HEADER ")
}

/// Convert an input file to uppercase FASTA, writing to the provided writer.
/// This is equivalent to running `any2fasta -q -u <file>` but in-process.
pub fn to_fasta_writer<W: Write>(path: &str, out: &mut W) -> Result<()> {
    let config = ConvConfig::new();

    let reader = open_input(path)
        .with_context(|| format!("Cannot open input file '{}'", path))?;

    let mut buf_reader = BufReader::new(reader);
    let mut first_line = String::new();
    if buf_reader.read_line(&mut first_line)? == 0 {
        anyhow::bail!("The input '{}' appears to be empty", path);
    }

    // Strip trailing newline
    if first_line.ends_with('\n') {
        first_line.pop();
        if first_line.ends_with('\r') {
            first_line.pop();
        }
    }

    let mut lines = Vec::new();
    lines.push(first_line);
    for line in buf_reader.lines() {
        lines.push(line?);
    }

    // Try each format detector
    let formats: Vec<Format<W>> = vec![
        Format { name: "GENBANK", detect: detect_genbank, parser: parse_genbank },
        Format { name: "UNIPROT", detect: detect_uniprot, parser: parse_embl },
        Format { name: "EMBL", detect: detect_embl, parser: parse_embl },
        Format { name: "FASTA", detect: detect_fasta, parser: parse_fasta },
        Format { name: "FASTQ", detect: detect_fastq, parser: parse_fastq },
        Format { name: "GFF", detect: detect_gff, parser: parse_gff },
        Format { name: "CLUSTAL", detect: detect_clustal, parser: parse_clustal },
        Format { name: "STOCKHOLM", detect: detect_stockholm, parser: parse_stockholm },
        Format { name: "GFA", detect: detect_gfa, parser: parse_gfa },
        Format { name: "PDB", detect: detect_pdb, parser: parse_pdb },
    ];

    for fmt in &formats {
        if (fmt.detect)(&lines[0]) {
            let nseq = (fmt.parser)(&config, &mut lines, out)?;
            if nseq == 0 {
                anyhow::bail!("No sequences found in '{}'", path);
            }
            return Ok(());
        }
    }

    anyhow::bail!(
        "Unfamiliar format with first line: {}",
        lines[0].trim()
    );
}

/// Convert an input file to uppercase FASTA string (in-memory).
pub fn to_fasta_string(path: &str) -> Result<String> {
    let mut buf = Vec::new();
    to_fasta_writer(path, &mut buf)?;
    Ok(String::from_utf8(buf)?)
}
