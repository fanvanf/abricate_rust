//! Database management — listing, setup, and metadata for abricate databases.
//! Each database is a subdirectory of the data directory containing a `sequences` FASTA file.

use crate::blast::{self, BlastPaths};
use anyhow::{Context, Result};
use std::fs;
use std::path::{Path, PathBuf};

/// Information about a single database
#[derive(Debug, Clone)]
pub struct DatabaseInfo {
    pub name: String,
    pub path: PathBuf,
    pub has_sequences: bool,
    pub has_blast_index: bool,
    pub db_type: String,
    pub num_sequences: usize,
}

/// Resolve the database data directory.
/// Priority: --datadir argument > ABRICATE_DATADIR env > default relative path
pub fn resolve_datadir(datadir_arg: Option<&str>) -> Result<PathBuf> {
    if let Some(d) = datadir_arg {
        let p = PathBuf::from(d);
        if !p.exists() {
            anyhow::bail!("Data directory does not exist: {}", p.display());
        }
        return Ok(p);
    }

    if let Ok(env_dir) = std::env::var("ABRICATE_DATADIR") {
        let p = PathBuf::from(&env_dir);
        if p.exists() {
            return Ok(p);
        }
    }

    // Default: try common locations relative to executable
    let candidates: Vec<Option<PathBuf>> = vec![
        // Relative to current directory
        Some(PathBuf::from("db")),
        // Relative to executable
        std::env::current_exe()
            .ok()
            .and_then(|e| e.parent().map(|p| p.join("db"))),
        std::env::current_exe()
            .ok()
            .and_then(|e| e.parent().and_then(|p| p.parent().map(|pp| pp.join("db")))),
        // User home
        std::env::var("HOME")
            .ok()
            .map(|h| Path::new(&h).join(".abricate").join("db")),
        std::env::var("USERPROFILE")
            .ok()
            .map(|h| Path::new(&h).join(".abricate").join("db")),
    ];

    for candidate in candidates.iter().flatten() {
        if candidate.exists() {
            return Ok(candidate.clone());
        }
    }

    // Last resort: return the default path (will fail on use)
    Ok(PathBuf::from("db"))
}

/// List all databases in the data directory
pub fn list_databases(datadir: &Path) -> Result<Vec<DatabaseInfo>> {
    let mut dbs = Vec::new();

    if !datadir.exists() {
        return Ok(dbs);
    }

    let entries = fs::read_dir(datadir)
        .with_context(|| format!("Cannot read data directory: {}", datadir.display()))?;

    for entry in entries {
        let entry = entry?;
        let path = entry.path();
        if !path.is_dir() {
            continue;
        }

        let name = entry.file_name().to_string_lossy().to_string();
        if name.starts_with('.') {
            continue;
        }

        let seq_file = path.join("sequences");
        let has_sequences = seq_file.exists();

        // Detect BLAST index files
        let has_blast_index = has_blast_index(&path, &seq_file);

        // Determine db type
        let db_type = if has_sequences {
            blast::detect_db_type(&seq_file).unwrap_or("nucl").to_string()
        } else {
            "unknown".to_string()
        };

        // Count sequences
        let num_sequences = if has_sequences {
            count_sequences(&seq_file).unwrap_or(0)
        } else {
            0
        };

        dbs.push(DatabaseInfo {
            name,
            path,
            has_sequences,
            has_blast_index,
            db_type,
            num_sequences,
        });
    }

    dbs.sort_by(|a, b| a.name.cmp(&b.name));
    Ok(dbs)
}

/// Check if BLAST index files exist for a database
fn has_blast_index(db_dir: &Path, seq_file: &Path) -> bool {
    // BLAST creates index files with the same prefix as the input file
    // For nucleotide: .nhr, .nin, .nsq
    // For protein: .phr, .pin, .psq
    let stem = seq_file.file_stem().map(|s| s.to_string_lossy().to_string());

    if let Some(stem) = stem {
        let nhr = db_dir.join(format!("{}.nhr", stem));
        let _nin = db_dir.join(format!("{}.nin", stem));
        let nsq = db_dir.join(format!("{}.nsq", stem));
        let phr = db_dir.join(format!("{}.phr", stem));
        let _pin = db_dir.join(format!("{}.pin", stem));
        let psq = db_dir.join(format!("{}.psq", stem));

        // Check for nucleotide index
        if nhr.exists() && nsq.exists() {
            return true;
        }
        // Check for protein index
        if phr.exists() && psq.exists() {
            return true;
        }
    }

    // Also check with the full sequences filename as prefix
    let nhr = db_dir.join("sequences.nhr");
    let nsq = db_dir.join("sequences.nsq");
    let phr = db_dir.join("sequences.phr");
    let psq = db_dir.join("sequences.psq");

    (nhr.exists() && nsq.exists()) || (phr.exists() && psq.exists())
}

/// Count the number of sequences in a FASTA file
fn count_sequences(path: &Path) -> Result<usize> {
    let content = fs::read_to_string(path)?;
    let count = content.lines().filter(|l| l.starts_with('>')).count();
    Ok(count)
}

/// Setup (build BLAST index) for all databases in the data directory
pub fn setup_all_databases(datadir: &Path, paths: &BlastPaths) -> Result<()> {
    let dbs = list_databases(datadir)?;

    if dbs.is_empty() {
        eprintln!("No databases found in {}", datadir.display());
        return Ok(());
    }

    for db in &dbs {
        if !db.has_sequences {
            eprintln!("[{}] SKIP - no sequences file", db.name);
            continue;
        }

        eprintln!("[{}] Building BLAST index (type: {})...", db.name, db.db_type);

        let seq_file = db.path.join("sequences");
        blast::make_blast_db(paths, &seq_file, &db.db_type)?;

        eprintln!("[{}] Done - {} sequences indexed", db.name, db.num_sequences);
    }

    eprintln!("All databases set up successfully.");
    Ok(())
}

/// Setup a single database
pub fn setup_database(datadir: &Path, db_name: &str, paths: &BlastPaths) -> Result<()> {
    let db_dir = datadir.join(db_name);
    if !db_dir.exists() {
        anyhow::bail!("Database '{}' not found in {}", db_name, datadir.display());
    }

    let seq_file = db_dir.join("sequences");
    if !seq_file.exists() {
        anyhow::bail!(
            "Database '{}' has no 'sequences' file in {}",
            db_name,
            db_dir.display()
        );
    }

    let db_type = blast::detect_db_type(&seq_file)?;
    eprintln!("[{}] Building BLAST index (type: {})...", db_name, db_type);

    blast::make_blast_db(paths, &seq_file, db_type)?;

    eprintln!("[{}] Done", db_name);
    Ok(())
}

/// Get database path for a specific database name
pub fn get_db_path(datadir: &Path, db_name: &str) -> Result<PathBuf> {
    let db_dir = datadir.join(db_name);
    if !db_dir.exists() {
        // List available databases for helpful error message
        let dbs = list_databases(datadir)?;
        let available: Vec<String> = dbs.iter().map(|d| d.name.clone()).collect();
        anyhow::bail!(
            "Database '{}' not found in {}. Available: {}",
            db_name,
            datadir.display(),
            if available.is_empty() {
                "(none)".to_string()
            } else {
                available.join(", ")
            }
        );
    }

    let seq_file = db_dir.join("sequences");
    if !seq_file.exists() {
        anyhow::bail!(
            "Database '{}' exists but has no 'sequences' file. Run --setupdb first.",
            db_name
        );
    }

    Ok(db_dir)
}
