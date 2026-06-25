//! BLAST+ binary management, database creation, and result parsing.
//! Handles cross-platform BLAST invocation with configurable binary paths.

use anyhow::{Context, Result};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// Separator used in abricate FASTA headers: `DB~~~GENE~~~ACCESSION~~~PRODUCT`
pub const ID_SEP: &str = "~~~";

/// BLAST tabular output fields used by abricate
const BLAST_OUTFMT: &str =
    "6 qseqid qstart qend qlen sseqid sstart send slen sstrand evalue length pident gaps gapopen stitle";

/// A single BLAST hit parsed from tabular output
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct BlastHit {
    pub qseqid: String,
    pub qstart: i64,
    pub qend: i64,
    pub qlen: i64,
    pub sseqid: String,
    pub sstart: i64,
    pub send: i64,
    pub slen: i64,
    pub sstrand: String,
    pub evalue: f64,
    pub length: i64,
    pub pident: f64,
    pub gaps: i64,
    pub gapopen: i64,
    pub stitle: String,
}

/// Manages BLAST+ binary paths across platforms
pub struct BlastPaths {
    pub blastn: PathBuf,
    pub blastx: PathBuf,
    pub makeblastdb: PathBuf,
    pub blastdbcmd: PathBuf,
}

impl BlastPaths {
    /// Resolve BLAST+ binary paths from --blastdir, env var, or PATH
    pub fn resolve(blastdir: Option<&str>) -> Result<Self> {
        let dir = blastdir
            .map(|s| s.to_string())
            .or_else(|| std::env::var("ABRICATE_BLAST_DIR").ok())
            .filter(|s| !s.is_empty());

        let exe_suffix = if cfg!(windows) { ".exe" } else { "" };

        let (blastn, blastx, makeblastdb, blastdbcmd) = if let Some(ref dir) = dir {
            let d = Path::new(dir);
            (
                d.join(format!("blastn{}", exe_suffix)),
                d.join(format!("blastx{}", exe_suffix)),
                d.join(format!("makeblastdb{}", exe_suffix)),
                d.join(format!("blastdbcmd{}", exe_suffix)),
            )
        } else {
            // Try to find in PATH
            (
                PathBuf::from(format!("blastn{}", exe_suffix)),
                PathBuf::from(format!("blastx{}", exe_suffix)),
                PathBuf::from(format!("makeblastdb{}", exe_suffix)),
                PathBuf::from(format!("blastdbcmd{}", exe_suffix)),
            )
        };

        // Verify at least blastn and makeblastdb are accessible
        let paths = BlastPaths {
            blastn,
            blastx,
            makeblastdb,
            blastdbcmd,
        };

        Ok(paths)
    }

    /// Check if all required BLAST binaries are accessible
    pub fn check(&self) -> Vec<(String, bool)> {
        let tools = [
            ("blastn", &self.blastn),
            ("blastx", &self.blastx),
            ("makeblastdb", &self.makeblastdb),
            ("blastdbcmd", &self.blastdbcmd),
        ];

        tools
            .iter()
            .map(|(name, path)| {
                let found = if path.is_absolute() {
                    path.exists()
                } else {
                    // Check if command is in PATH
                    Command::new(path)
                        .arg("-help")
                        .stdout(Stdio::null())
                        .stderr(Stdio::null())
                        .status()
                        .is_ok()
                };
                (name.to_string(), found)
            })
            .collect()
    }

    /// Get a display string for the binary directory
    pub fn dir_display(&self) -> String {
        if self.blastn.is_absolute() {
            self.blastn
                .parent()
                .map(|p| p.display().to_string())
                .unwrap_or_else(|| "(unknown)".to_string())
        } else {
            "PATH".to_string()
        }
    }
}

/// Determine if a sequences file is nucleotide or protein by examining content
pub fn detect_db_type(sequences_path: &Path) -> Result<&'static str> {
    let content = std::fs::read_to_string(sequences_path)
        .with_context(|| format!("Cannot read sequences file: {}", sequences_path.display()))?;

    let mut total_len = 0usize;
    let mut non_atgc = 0usize;

    for line in content.lines() {
        if line.starts_with('>') {
            continue;
        }
        let trimmed = line.trim();
        for c in trimmed.chars() {
            if !c.is_whitespace() {
                total_len += 1;
                if !"ATCGNatcgn".contains(c) {
                    non_atgc += 1;
                }
            }
        }
    }

    if total_len == 0 {
        anyhow::bail!("No sequence data found in {}", sequences_path.display());
    }

    if non_atgc as f64 / total_len as f64 > 0.5 {
        Ok("prot")
    } else {
        Ok("nucl")
    }
}

/// Build a BLAST database from a sequences file using makeblastdb
pub fn make_blast_db(paths: &BlastPaths, sequences_path: &Path, db_type: &str) -> Result<()> {
    let dbtype = if db_type == "prot" { "prot" } else { "nucl" };

    // makeblastdb creates files with the same prefix as the input (without extension)
    // We use the sequences path itself as the db prefix
    let output = Command::new(&paths.makeblastdb)
        .arg("-in")
        .arg(sequences_path)
        .arg("-dbtype")
        .arg(dbtype)
        .arg("-parse_seqids")
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .with_context(|| {
            format!(
                "Failed to execute makeblastdb at {}",
                paths.makeblastdb.display()
            )
        })?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!("makeblastdb failed: {}", stderr);
    }

    Ok(())
}

/// Get database type using blastdbcmd -info, falling back to content detection
pub fn get_db_type(paths: &BlastPaths, db_path: &Path) -> Result<String> {
    // Try blastdbcmd -info first
    let output = Command::new(&paths.blastdbcmd)
        .arg("-db")
        .arg(db_path)
        .arg("-info")
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output();

    if let Ok(out) = output {
        if out.status.success() {
            let info = String::from_utf8_lossy(&out.stdout);
            // blastdbcmd -info output contains "Database: ... (type: Nucleotide)" or "(type: Protein)"
            if info.contains("Nucleotide") {
                return Ok("nucl".to_string());
            } else if info.contains("Protein") {
                return Ok("prot".to_string());
            }
        }
    }

    // Fallback: detect from sequences file content
    let sequences = db_path.join("sequences");
    if sequences.exists() {
        let detected = detect_db_type(&sequences)?;
        Ok(detected.to_string())
    } else {
        // Default to nucleotide
        Ok("nucl".to_string())
    }
}

/// Run BLAST (blastn or blastx) with the given FASTA input piped via stdin
pub fn run_blast(
    paths: &BlastPaths,
    db_path: &Path,
    db_type: &str,
    fasta_input: &str,
    minid: f64,
    threads: usize,
) -> Result<Vec<BlastHit>> {
    let blast_cmd = if db_type == "nucl" {
        &paths.blastn
    } else {
        &paths.blastx
    };

    let mut cmd = Command::new(blast_cmd);
    cmd.arg("-db").arg(db_path.join("sequences"));

    if db_type == "nucl" {
        cmd.arg("-task").arg("blastn");
        cmd.arg("-dust").arg("no");
        cmd.arg("-perc_identity").arg(format!("{}", minid as u64));
    } else {
        cmd.arg("-task").arg("blastx-fast");
        cmd.arg("-seg").arg("no");
    }

    cmd.arg("-outfmt").arg(BLAST_OUTFMT);
    cmd.arg("-num_threads").arg(format!("{}", threads));
    cmd.arg("-evalue").arg("1E-20");
    cmd.arg("-culling_limit").arg("5");
    cmd.arg("-max_target_seqs").arg("10000");

    cmd.stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());

    let mut child = cmd
        .spawn()
        .with_context(|| format!("Failed to spawn BLAST process: {}", blast_cmd.display()))?;

    // Write FASTA to stdin
    if let Some(mut stdin) = child.stdin.take() {
        stdin
            .write_all(fasta_input.as_bytes())
            .context("Failed to write to BLAST stdin")?;
        // stdin is dropped here, closing the pipe
    }

    let output = child
        .wait_with_output()
        .context("Failed to wait for BLAST process")?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        // BLAST sometimes prints warnings to stderr even on success
        if stderr.contains("Error") || stderr.contains("ERROR") {
            anyhow::bail!("BLAST failed: {}", stderr);
        }
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    parse_blast_output(&stdout)
}

/// Parse BLAST tabular output into BlastHit structs
fn parse_blast_output(output: &str) -> Result<Vec<BlastHit>> {
    let mut hits = Vec::new();

    for line in output.lines() {
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 15 {
            continue;
        }

        let hit = BlastHit {
            qseqid: fields[0].to_string(),
            qstart: fields[1].parse().unwrap_or(0),
            qend: fields[2].parse().unwrap_or(0),
            qlen: fields[3].parse().unwrap_or(0),
            sseqid: fields[4].to_string(),
            sstart: fields[5].parse().unwrap_or(0),
            send: fields[6].parse().unwrap_or(0),
            slen: fields[7].parse().unwrap_or(0),
            sstrand: fields[8].to_string(),
            evalue: fields[9].parse().unwrap_or(0.0),
            length: fields[10].parse().unwrap_or(0),
            pident: fields[11].parse().unwrap_or(0.0),
            gaps: fields[12].parse().unwrap_or(0),
            gapopen: fields[13].parse().unwrap_or(0),
            stitle: fields[14].to_string(),
        };

        hits.push(hit);
    }

    Ok(hits)
}

/// Generate an ASCII coverage map for a hit
/// Uses '=' for covered positions and '.' for uncovered positions on the subject sequence
pub fn coverage_map(sstart: i64, send: i64, slen: i64) -> String {
    if slen <= 0 {
        return String::new();
    }
    let start = sstart.min(send);
    let end = sstart.max(send);
    let mut map = String::with_capacity(slen as usize);
    for i in 1..=slen {
        if i >= start && i <= end {
            map.push('=');
        } else {
            map.push('.');
        }
    }
    map
}

/// Parse a subject sequence ID (sseqid or stitle) into (database, gene, accession, product)
/// Format: `DB~~~GENE~~~ACCESSION~~~PRODUCT`
pub fn parse_sseqid(stitle: &str) -> (String, String, String, String) {
    let parts: Vec<&str> = stitle.splitn(4, ID_SEP).collect();
    match parts.len() {
        4 => (
            parts[0].to_string(),
            parts[1].to_string(),
            parts[2].to_string(),
            parts[3].to_string(),
        ),
        3 => (
            parts[0].to_string(),
            parts[1].to_string(),
            parts[2].to_string(),
            String::new(),
        ),
        2 => (
            parts[0].to_string(),
            parts[1].to_string(),
            String::new(),
            String::new(),
        ),
        1 => (parts[0].to_string(), String::new(), String::new(), String::new()),
        _ => (String::new(), String::new(), String::new(), String::new()),
    }
}
