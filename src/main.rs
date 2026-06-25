//! abricate_rust — Rust implementation of ABRicate
//! Mass screening of contigs for antimicrobial resistance or virulence genes.
//! Cross-platform (Windows, Linux, macOS) with integrated any2fasta and configurable BLAST+.

mod blast;
mod database;
mod fasta;
mod getdb;
mod report;
mod summary;

use anyhow::{Context, Result};
use clap::Parser;

// ---------------------------------------------------------------------------
// CLI definition
// ---------------------------------------------------------------------------
#[derive(Parser, Debug)]
#[command(
    name = "abricate",
    version,
    about = "Mass screening of contigs for antimicrobial resistance or virulence genes (Rust implementation)"
)]
struct Args {
    /// List all available databases in the data directory
    #[arg(long)]
    list: bool,

    /// Build BLAST indices for all databases (or a specific one with --db)
    #[arg(long)]
    setupdb: bool,

    /// Check that all required dependencies are available
    #[arg(long)]
    check: bool,

    /// Generate a summary matrix from previously generated abricate reports
    #[arg(long)]
    summary: bool,

    /// Use %IDENTITY instead of %COVERAGE in summary mode
    #[arg(long)]
    identity: bool,

    /// Download and set up a database from the internet.
    /// Available: resfinder, plasmidfinder, megares, argannot, card, ncbi,
    /// vfdb, ecoli_vf, upec_expec_vf, ecoh, bacmet2, victors
    /// Use 'list' to show available databases
    #[arg(long, value_name = "DATABASE")]
    getdb: Option<String>,

    /// Force re-download even if database already exists
    #[arg(long)]
    force: bool,

    /// Database to use [default: ncbi]
    #[arg(long, default_value = "ncbi")]
    db: String,

    /// Data directory containing databases [default: auto-detect]
    #[arg(long)]
    datadir: Option<String>,

    /// Directory containing BLAST+ binaries (blastn, makeblastdb, etc.)
    /// Can also be set via ABRICATE_BLAST_DIR environment variable.
    /// On Linux/macOS: optional (will use PATH if not specified).
    /// On Windows: required if BLAST+ is not in PATH.
    #[arg(long)]
    blastdir: Option<String>,

    /// Minimum DNA %identity [default: 80]
    #[arg(long, default_value_t = 80.0)]
    minid: f64,

    /// Minimum DNA %coverage [default: 80]
    #[arg(long, default_value_t = 80.0)]
    mincov: f64,

    /// Number of BLAST threads [default: 1]
    #[arg(long, short = 't', default_value_t = 1)]
    threads: usize,

    /// Output format: tsv, csv, bed, gff, json [default: tsv]
    #[arg(long, default_value = "tsv")]
    outfmt: String,

    /// Suppress header row in output
    #[arg(long)]
    noheader: bool,

    /// Strip directory paths from FILE column
    #[arg(long)]
    nopath: bool,

    /// File containing a list of input files (one per line)
    #[arg(long)]
    fofn: Option<String>,

    /// Input files (FASTA/GenBank/EMBL/GFF/etc., optionally .gz/.bz2/.zip compressed)
    files: Vec<String>,
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------
fn main() {
    let args = Args::parse();

    if let Err(e) = run(args) {
        eprintln!("ERROR: {}", e);
        std::process::exit(1);
    }
}

fn run(args: Args) -> Result<()> {
    // Handle special modes first
    if args.check {
        return run_check(&args);
    }

    if args.list {
        return run_list(&args);
    }

    // --getdb: download and set up a database
    if let Some(db_name) = &args.getdb {
        if db_name == "list" {
            eprintln!("Available databases for download:");
            for db in getdb::AVAILABLE_DATABASES {
                eprintln!("  {}", db);
            }
            return Ok(());
        }
        return run_getdb(&args, db_name);
    }

    if args.setupdb {
        return run_setupdb(&args);
    }

    if args.summary {
        return run_summary(&args);
    }

    // Default: run detection pipeline
    run_detect(&args)
}

// ---------------------------------------------------------------------------
// --check: verify dependencies
// ---------------------------------------------------------------------------
fn run_check(args: &Args) -> Result<()> {
    eprintln!("Checking dependencies for abricate_rust v{}", env!("CARGO_PKG_VERSION"));
    eprintln!();

    // Check BLAST+
    let paths = blast::BlastPaths::resolve(args.blastdir.as_deref())?;
    eprintln!("BLAST+ binaries directory: {}", paths.dir_display());
    let checks = paths.check();
    let mut all_ok = true;
    for (name, found) in &checks {
        let status = if *found { "OK" } else { "MISSING" };
        eprintln!("  {:20} {}", name, status);
        if !found {
            all_ok = false;
        }
    }
    eprintln!();

    if !all_ok {
        eprintln!("WARNING: Some BLAST+ binaries are missing!");
        eprintln!("  Set --blastdir <PATH> or ABRICATE_BLAST_DIR environment variable");
        eprintln!("  to point to your BLAST+ installation directory.");
        eprintln!("  Example: abricate --blastdir /path/to/ncbi-blast/bin --check");
    }

    // Check data directory
    let datadir = database::resolve_datadir(args.datadir.as_deref())?;
    eprintln!("Data directory: {}", datadir.display());
    if datadir.exists() {
        let dbs = database::list_databases(&datadir)?;
        eprintln!("  Found {} databases: {}", dbs.len(), dbs.iter().map(|d| d.name.clone()).collect::<Vec<_>>().join(", "));
    } else {
        eprintln!("  WARNING: Data directory does not exist!");
    }
    eprintln!();

    // any2fasta is integrated, always available
    eprintln!("  {:20} OK (integrated)", "any2fasta");

    eprintln!();
    if all_ok {
        eprintln!("All dependencies satisfied.");
    } else {
        eprintln!("Some dependencies are missing. Please fix before proceeding.");
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// --list: list databases
// ---------------------------------------------------------------------------
fn run_list(args: &Args) -> Result<()> {
    let datadir = database::resolve_datadir(args.datadir.as_deref())?;

    if !datadir.exists() {
        eprintln!("Data directory not found: {}", datadir.display());
        eprintln!("Set --datadir to point to your database directory.");
        return Ok(());
    }

    let dbs = database::list_databases(&datadir)?;

    if dbs.is_empty() {
        eprintln!("No databases found in {}", datadir.display());
        return Ok(());
    }

    println!("{:<20} {:<6} {:>10}  {}", "DATABASE", "TYPE", "SEQUENCES", "BLAST_INDEX");
    println!("{}", "-".repeat(60));
    for db in &dbs {
        let idx_status = if db.has_blast_index { "YES" } else { "NO" };
        println!(
            "{:<20} {:<6} {:>10}  {}",
            db.name, db.db_type, db.num_sequences, idx_status
        );
    }
    println!();
    eprintln!("{} database(s) found.", dbs.len());

    Ok(())
}

// ---------------------------------------------------------------------------
// --setupdb: build BLAST indices
// ---------------------------------------------------------------------------
fn run_setupdb(args: &Args) -> Result<()> {
    let datadir = database::resolve_datadir(args.datadir.as_deref())?;
    let paths = blast::BlastPaths::resolve(args.blastdir.as_deref())?;
    paths.require_blast()?;

    if args.db != "ncbi" {
        // Setup specific database
        database::setup_database(&datadir, &args.db, &paths)?;
    } else {
        // Setup all databases
        database::setup_all_databases(&datadir, &paths)?;
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// --getdb: download and set up a database
// ---------------------------------------------------------------------------
fn run_getdb(args: &Args, db_name: &str) -> Result<()> {
    // For getdb, create the datadir if it doesn't exist
    let datadir = if let Some(d) = &args.datadir {
        let p = std::path::PathBuf::from(d);
        if !p.exists() {
            std::fs::create_dir_all(&p)?;
        }
        p
    } else if let Ok(env_dir) = std::env::var("ABRICATE_DATADIR") {
        let p = std::path::PathBuf::from(&env_dir);
        if !p.exists() {
            std::fs::create_dir_all(&p)?;
        }
        p
    } else {
        let p = database::resolve_datadir(None)?;
        if !p.exists() {
            std::fs::create_dir_all(&p)?;
        }
        p
    };

    let paths = blast::BlastPaths::resolve(args.blastdir.as_deref())?;
    paths.require_blast()?;

    getdb::download_database(db_name, &datadir, args.force, &paths)?;

    Ok(())
}

// ---------------------------------------------------------------------------
// --summary: generate summary matrix
// ---------------------------------------------------------------------------
fn run_summary(args: &Args) -> Result<()> {
    let files = collect_input_files(&args)?;
    if files.is_empty() {
        anyhow::bail!("No input files specified. Provide abricate TSV report files or use --fofn.");
    }

    let mut stdout = std::io::stdout().lock();
    summary::generate_summary(&files, args.identity, args.noheader, &mut stdout)?;

    Ok(())
}

// ---------------------------------------------------------------------------
// Default mode: run detection pipeline
// ---------------------------------------------------------------------------
fn run_detect(args: &Args) -> Result<()> {
    let files = collect_input_files(&args)?;
    if files.is_empty() {
        anyhow::bail!(
            "No input files specified.\n\
             Usage: abricate [OPTIONS] <FASTA_FILES...>\n\
             Use --fofn to specify a file containing input file paths."
        );
    }

    // Resolve paths
    let datadir = database::resolve_datadir(args.datadir.as_deref())?;
    let paths = blast::BlastPaths::resolve(args.blastdir.as_deref())?;
    paths.require_blast()?;

    // Get database path
    let db_path = database::get_db_path(&datadir, &args.db)?;

    // Determine database type
    let db_type = blast::get_db_type(&paths, &db_path)?;
    eprintln!("Using database: {} (type: {})", args.db, db_type);
    eprintln!("Parameters: minid={}, mincov={}, threads={}", args.minid, args.mincov, args.threads);
    eprintln!("Processing {} file(s)...", files.len());

    // Parse output format
    let format = report::OutFormat::parse(&args.outfmt)
        .map_err(|e| anyhow::anyhow!(e))?;

    let mut stdout = std::io::stdout().lock();
    let mut header_written = false;

    for (i, file) in files.iter().enumerate() {
        eprintln!("[{}/{}] Processing: {}", i + 1, files.len(), file);

        // Convert input to FASTA (integrated any2fasta, uppercase)
        let fasta_input = match fasta::to_fasta_string(file) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("  WARNING: Skipping '{}': {}", file, e);
                continue;
            }
        };

        if fasta_input.trim().is_empty() {
            eprintln!("  WARNING: No sequences in '{}'", file);
            continue;
        }

        // Run BLAST
        let hits = match blast::run_blast(
            &paths,
            &db_path,
            &db_type,
            &fasta_input,
            args.minid,
            args.threads,
        ) {
            Ok(h) => h,
            Err(e) => {
                eprintln!("  WARNING: BLAST failed for '{}': {}", file, e);
                continue;
            }
        };

        // Process hits into report rows
        let rows = report::process_hits(
            hits,
            file,
            &args.db,
            args.minid,
            args.mincov,
            args.nopath,
        );

        eprintln!("  Found {} hit(s)", rows.len());

        // Write output (header only once for all files)
        let suppress_header = args.noheader || header_written;
        report::write_report(&rows, format, suppress_header, &mut stdout)?;
        if !suppress_header {
            header_written = true;
        }
    }

    eprintln!("Done.");
    Ok(())
}

// ---------------------------------------------------------------------------
// Helper: collect input files from args and --fofn
// ---------------------------------------------------------------------------
fn collect_input_files(args: &Args) -> Result<Vec<String>> {
    let mut files = args.files.clone();

    if let Some(fofn) = &args.fofn {
        let content = std::fs::read_to_string(fofn)
            .with_context(|| format!("Cannot read FOFN file: {}", fofn))?;
        for line in content.lines() {
            let line = line.trim();
            if !line.is_empty() && !line.starts_with('#') {
                files.push(line.to_string());
            }
        }
    }

    Ok(files)
}
