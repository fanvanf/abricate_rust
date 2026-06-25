//! Database download module — integrated abricate-get_db functionality.
//! Downloads and parses AMR/virulence databases from various sources,
//! formats them as abricate FASTA (with ~~~ headers), and builds BLAST indices.

use crate::blast::{self, BlastPaths};
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs;
use std::io::{Read, Write};
use std::path::Path;
use std::process::Command;
use std::time::Duration;

/// All databases available for download
pub const AVAILABLE_DATABASES: &[&str] = &[
    "resfinder",
    "plasmidfinder",
    "megares",
    "argannot",
    "card",
    "ncbi",
    "vfdb",
    "ecoli_vf",
    "upec_expec_vf",
    "ecoh",
    "bacmet2",
    "victors",
];

const ABX_SEP: &str = ";";

// ===========================================================================
// Data structures
// ===========================================================================

/// A single sequence record from a database
#[derive(Debug, Clone)]
struct DbSeq {
    id: String,
    acc: String,
    desc: String,
    seq: String,
    abx: Vec<String>,
}

/// A raw FASTA record
struct FastaRecord {
    id: String,
    desc: String,
    seq: String,
}

// ===========================================================================
// Utility functions
// ===========================================================================

/// Download a URL to a file, skipping if already exists (unless force)
fn download(url: &str, dest: &Path, force: bool) -> Result<()> {
    if dest.exists() && !force {
        eprintln!("  Won't re-download existing {} (use --force)", dest.display());
        return Ok(());
    }

    eprintln!("  Downloading: {}", url);
    let resp = ureq::get(url)
        .timeout(Duration::from_secs(600))
        .call()
        .map_err(|e| anyhow::anyhow!("Download failed: {}", e))?;

    let mut file = fs::File::create(dest)
        .with_context(|| format!("Cannot create file: {}", dest.display()))?;
    let mut reader = resp.into_reader();
    std::io::copy(&mut reader, &mut file)?;

    let size = fs::metadata(dest)?.len();
    eprintln!("  Saved: {} ({} bytes)", dest.display(), size);
    Ok(())
}

/// Run git clone
fn run_git_clone(url: &str, dest: &Path, force: bool) -> Result<()> {
    if dest.exists() && !force {
        eprintln!("  Won't overwrite existing {} (use --force)", dest.display());
        return Ok(());
    }

    if dest.exists() {
        eprintln!("  Removing existing {}", dest.display());
        fs::remove_dir_all(dest)?;
    }

    eprintln!("  Cloning: {}", url);
    let output = Command::new("git")
        .arg("clone")
        .arg("--quiet")
        .arg(url)
        .arg(dest)
        .output()
        .context("Failed to run git clone. Is git installed?")?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!("git clone failed: {}", stderr);
    }

    Ok(())
}

/// Parse a FASTA string into records
fn parse_fasta(content: &str) -> Vec<FastaRecord> {
    let mut records = Vec::new();
    let mut current_id = String::new();
    let mut current_desc = String::new();
    let mut current_seq = String::new();
    let mut in_seq = false;

    for line in content.lines() {
        let line = line.trim_end();
        if line.starts_with('>') {
            if in_seq && !current_seq.is_empty() {
                records.push(FastaRecord {
                    id: current_id.clone(),
                    desc: current_desc.clone(),
                    seq: current_seq.clone(),
                });
            }
            let header = &line[1..];
            let mut parts = header.splitn(2, |c: char| c == ' ' || c == '\t');
            current_id = parts.next().unwrap_or("").to_string();
            current_desc = parts.next().unwrap_or("").trim().to_string();
            current_seq.clear();
            in_seq = true;
        } else if in_seq {
            let trimmed = line.trim();
            if !trimmed.is_empty() {
                current_seq.push_str(trimmed);
            }
        }
    }

    if in_seq && !current_seq.is_empty() {
        records.push(FastaRecord {
            id: current_id,
            desc: current_desc,
            seq: current_seq,
        });
    }

    records
}

/// Load a FASTA file into DbSeq records (with sanitization)
fn load_fasta_file(path: &Path) -> Result<Vec<DbSeq>> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("Cannot read FASTA file: {}", path.display()))?;
    let raw = parse_fasta(&content);

    let mut seqs = Vec::new();
    let mut seen: HashMap<String, usize> = HashMap::new();

    for rec in raw {
        if rec.id.is_empty() {
            continue;
        }

        // Handle duplicate IDs
        let id = if let Some(_count) = seen.get(&rec.id) {
            let new_id = format!("{}_dupe", rec.id);
            eprintln!("  WARNING: Duplicate ID '{}' in {}", rec.id, path.display());
            new_id
        } else {
            rec.id.clone()
        };
        seen.insert(rec.id, 1);

        // Sanitize sequence: uppercase, replace invalid chars
        let s = rec.seq.to_uppercase();
        let is_dna = is_dna_seq(&s);
        let sanitized = if is_dna {
            s.chars()
                .map(|c| if "AGTCN".contains(c) { c } else { 'N' })
                .collect()
        } else {
            s.chars()
                .map(|c| if c.is_ascii_uppercase() { c } else { 'X' })
                .collect()
        };

        seqs.push(DbSeq {
            id,
            acc: String::new(),
            desc: rec.desc,
            seq: sanitized,
            abx: Vec::new(),
        });
    }

    eprintln!("  load_fasta: read {} sequences from {}", seqs.len(), path.display());
    Ok(seqs)
}

/// Check if a sequence is DNA
fn is_dna_seq(seq: &str) -> bool {
    if seq.is_empty() {
        return true;
    }
    let non_atgc = seq.chars().filter(|c| !"ATCGN".contains(*c)).count();
    (non_atgc as f64) / (seq.len() as f64) < 0.3
}

/// Save sequences in abricate FASTA format
fn save_fasta(path: &Path, db_name: &str, seqs: &[DbSeq]) -> Result<()> {
    eprintln!("  save_fasta: writing {} sequences to {}", seqs.len(), path.display());
    let mut output = fs::File::create(path)?;

    for s in seqs {
        let abx_str = s
            .abx
            .iter()
            .map(|a| a.replace(' ', "_"))
            .collect::<Vec<_>>()
            .join(ABX_SEP);
        let header_id = vec![db_name, &s.id, &s.acc, &abx_str].join("~~~");
        let desc = if s.desc.is_empty() { &s.id } else { &s.desc };
        writeln!(output, ">{} {}", header_id, desc)?;
        writeln!(output, "{}", s.seq)?;
    }

    Ok(())
}

/// Remove duplicate sequences (by sequence content)
fn dedupe_seq(seqs: Vec<DbSeq>) -> Vec<DbSeq> {
    let total = seqs.len();
    let mut seen: HashMap<String, String> = HashMap::new();
    let mut good = Vec::new();

    for s in seqs {
        if let Some(existing) = seen.get(&s.seq) {
            eprintln!(
                "  WARNING: duplicate {}bp sequence: {} ~ {}",
                s.seq.len(),
                s.id,
                existing
            );
        } else {
            good.push(s.clone());
        }
        seen.insert(s.seq, s.id);
    }

    eprintln!("  dedupe: read {} / kept {}", total, good.len());
    good
}

/// Load a TSV file into a HashMap keyed by a specific column
fn load_tabular(path: &Path, key_col: usize) -> Result<HashMap<String, HashMap<String, String>>> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("Cannot read TSV file: {}", path.display()))?;
    let mut result = HashMap::new();
    let mut headers: Vec<String> = Vec::new();

    for (line_num, line) in content.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();

        if headers.is_empty() {
            headers = cols.iter().map(|s| s.to_string()).collect();
            continue;
        }

        if cols.len() < headers.len() {
            eprintln!(
                "  WARNING: row {} has {} cols, expected {}",
                line_num + 1,
                cols.len(),
                headers.len()
            );
            continue;
        }

        let key = cols.get(key_col).unwrap_or(&"").to_string();
        if key.is_empty() {
            continue;
        }

        let row: HashMap<String, String> = headers
            .iter()
            .zip(cols.iter())
            .map(|(h, c)| (h.clone(), c.to_string()))
            .collect();

        result.entry(key).or_insert(row);
    }

    Ok(result)
}

/// Extract all files from a ZIP archive
fn extract_zip(zip_path: &Path, dest_dir: &Path) -> Result<()> {
    let file = fs::File::open(zip_path)?;
    let mut archive = zip::ZipArchive::new(file)?;

    for i in 0..archive.len() {
        let mut entry = archive.by_index(i)?;
        let name = entry.name().to_string();
        if name.ends_with('/') {
            continue;
        }

        let filename = Path::new(&name)
            .file_name()
            .ok_or_else(|| anyhow::anyhow!("Invalid filename in zip: {}", name))?;
        let outpath = dest_dir.join(filename);

        let mut out = fs::File::create(&outpath)?;
        std::io::copy(&mut entry, &mut out)?;
    }

    Ok(())
}

/// Extract a specific file from a tar.bz2 archive
fn extract_file_from_tar_bz2(tar_path: &Path, target_filename: &str) -> Result<Vec<u8>> {
    let file = fs::File::open(tar_path)?;
    let decoder = bzip2::read::BzDecoder::new(file);
    let mut archive = tar::Archive::new(decoder);

    for entry in archive.entries()? {
        let mut entry = entry?;
        let path = entry.path()?;
        let path_str = path.to_string_lossy();
        if path_str.ends_with(target_filename) {
            let mut buf = Vec::new();
            entry.read_to_end(&mut buf)?;
            return Ok(buf);
        }
    }

    anyhow::bail!("File '{}' not found in tar archive", target_filename)
}

/// Decompress a gzip file
fn decompress_gz(gz_path: &Path) -> Result<String> {
    let file = fs::File::open(gz_path)?;
    let mut decoder = flate2::read::GzDecoder::new(file);
    let mut content = String::new();
    decoder.read_to_string(&mut content)?;
    Ok(content)
}

// ===========================================================================
// Database-specific download and parse functions
// ===========================================================================

/// ResFinder database
fn get_resfinder(workdir: &Path, force: bool) -> Result<Vec<DbSeq>> {
    let name = "resfinder_db";
    let url = "https://bitbucket.org/genomicepidemiology/resfinder_db.git";
    let repo_dir = workdir.join(name);

    run_git_clone(url, &repo_dir, force)?;

    // Parse phenotypes.txt for antibiotic annotations
    let metafn = repo_dir.join("phenotypes.txt");
    let mut anno: HashMap<String, Vec<String>> = HashMap::new();

    if metafn.exists() {
        let content = fs::read_to_string(&metafn)?;
        for line in content.lines() {
            if line.starts_with('#') {
                continue;
            }
            let x: Vec<&str> = line.split('\t').collect();
            if x.is_empty() || x[0].is_empty() {
                continue;
            }
            // Parse gene name from first column: "aac(6')-Ia_1_X04555" -> "aac(6')-Ia"
            let gene = match x[0].rsplitn(2, '_').nth(1) {
                Some(g) => g.to_string(),
                None => x[0].to_string(),
            };

            // Column 2 (0-indexed) is the phenotype
            if x.len() > 2 {
                let abx: Vec<String> = x[2]
                    .split(", ")
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty() && !s.to_lowercase().contains("unknown") && !s.to_lowercase().contains("notes") && !s.to_lowercase().contains("none"))
                    .collect();
                if !abx.is_empty() {
                    anno.entry(gene).or_default().extend(abx);
                }
            }
        }
    }
    eprintln!("  resfinder: {} annotated genes", anno.len());

    // Parse .fsa files
    let mut all_seqs = Vec::new();
    for entry in fs::read_dir(&repo_dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map(|e| e == "fsa").unwrap_or(false) {
            // Fix broken FASTA: insert newline before '>' if preceded by a letter
            let content = fs::read_to_string(&path)?;
            let fixed = regex_replace_inline_gt(&content);
            let raw = parse_fasta(&fixed);

            for rec in raw {
                let s = rec.seq.to_uppercase();
                let sanitized: String = s.chars().map(|c| if "AGTCN".contains(c) { c } else { 'N' }).collect();

                // Parse ID: "aac(6')-Ib_2_M23634" -> id="aac(6')-Ib", copy=2, acc="M23634"
                let (gene_id, copy, acc) = parse_resfinder_id(&rec.id);

                let abx = anno.get(&gene_id).cloned().unwrap_or_default();

                all_seqs.push(DbSeq {
                    id: format!("{}_{}", gene_id, copy),
                    acc,
                    desc: gene_id,
                    seq: sanitized,
                    abx,
                });
            }
        }
    }

    Ok(all_seqs)
}

/// Fix inline '>' in FASTA (like Perl sed -i.bak 's/\([A-Z]\)>/\1\n>/gi')
fn regex_replace_inline_gt(content: &str) -> String {
    let mut result = String::with_capacity(content.len());
    let chars: Vec<char> = content.chars().collect();
    for i in 0..chars.len() {
        if chars[i] == '>' && i > 0 && chars[i - 1].is_ascii_alphabetic() {
            result.push('\n');
        }
        result.push(chars[i]);
    }
    result
}

/// Parse ResFinder ID: "aac(6')-Ib_2_M23634" -> ("aac(6')-Ib", "2", "M23634")
fn parse_resfinder_id(id: &str) -> (String, String, String) {
    let parts: Vec<&str> = id.rsplitn(3, '_').collect();
    match parts.len() {
        3 => (parts[2].to_string(), parts[1].to_string(), parts[0].to_string()),
        2 => (parts[1].to_string(), parts[0].to_string(), String::new()),
        _ => (id.to_string(), "1".to_string(), String::new()),
    }
}

/// PlasmidFinder database
fn get_plasmidfinder(workdir: &Path, force: bool) -> Result<Vec<DbSeq>> {
    let zip = workdir.join("plasmidfinder.zip");
    download(
        "https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/HEAD.zip",
        &zip,
        force,
    )?;
    extract_zip(&zip, workdir)?;

    let mut all_seqs = Vec::new();
    for entry in fs::read_dir(workdir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map(|e| e == "fsa").unwrap_or(false) {
            let seqs = load_fasta_file(&path)?;
            for mut s in seqs {
                // Parse ID: "IncFIA_1_AP001918.1" -> id="IncFIA", acc="AP001918.1"
                let desc = s.id.clone();
                let id_str = s.id.clone();
                if let Some(pos) = id_str.rfind('_') {
                    let acc_part = &id_str[pos + 1..];
                    let rest = &id_str[..pos];
                    if let Some(pos2) = rest.rfind('_') {
                        let id = rest[..pos2].trim_end_matches('_').to_string();
                        s.id = if id.is_empty() { desc.clone() } else { id };
                        s.acc = acc_part.to_string();
                    }
                }
                s.desc = desc;
                all_seqs.push(s);
            }
        }
    }

    Ok(all_seqs)
}

/// MEGARes database
fn get_megares(workdir: &Path, force: bool) -> Result<Vec<DbSeq>> {
    let zip = workdir.join("megares.zip");
    download("https://www.meglab.org/downloads/megares_v3.00.zip", &zip, force)?;
    extract_zip(&zip, workdir)?;

    // Find the drugs fasta file
    let mut fasta_path = None;
    for entry in fs::read_dir(workdir)? {
        let entry = entry?;
        let path = entry.path();
        let name = path.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default();
        if name.starts_with("megares_drugs_") && path.extension().map(|e| e == "fasta" || e == "fa").unwrap_or(false) {
            fasta_path = Some(path);
            break;
        }
    }

    let fasta_path = fasta_path.ok_or_else(|| anyhow::anyhow!("megares drugs fasta not found after extraction"))?;
    let seqs = load_fasta_file(&fasta_path)?;

    let mut okseq = Vec::new();
    for mut s in seqs {
        // Parse: >MEG_372|Multi-compound|Drug_and_biocide_resistance|...|ABEM
        let id_str = s.id.clone();
        let parts: Vec<&str> = id_str.split('|').collect();
        if parts.len() < 5 {
            continue;
        }
        let id = parts[0];
        let dbtype = parts[1];
        let class = parts[2];
        let mech = parts[3];
        let group = parts[4];
        let note = parts.get(5).map(|s| s.to_string()).unwrap_or_default();

        if !note.is_empty() {
            eprintln!("  Skipping {} due to: {}", id, note);
            continue;
        }

        s.id = group.to_string();
        s.acc = id.to_string();
        s.desc = format!("{}:{}:{}", dbtype, class, mech);
        okseq.push(s);
    }

    Ok(okseq)
}

/// ARG-ANNOT database
fn get_argannot(workdir: &Path, force: bool) -> Result<Vec<DbSeq>> {
    let fasta = workdir.join("arg-annot.fa");
    download(
        "https://www.mediterranee-infection.com/wp-content/uploads/2019/09/ARG-ANNOT_NT_V6_July2019.txt",
        &fasta,
        force,
    )?;

    // Fix backslashes in FASTA
    let content = fs::read_to_string(&fasta)?;
    let fixed = content.replace('\\', "");
    fs::write(&fasta, &fixed)?;

    let seqs = load_fasta_file(&fasta)?;

    let mut result = Vec::new();
    for mut s in seqs {
        // Parse: ">(AGly)Aac2-Ie:NC_011896:3039059-3039607:549"
        let id_str = s.id.clone();
        let parts: Vec<&str> = id_str.split(':').collect();
        if parts.len() >= 3 {
            s.id = parts[0].to_string();
            s.acc = format!("{}:{}", parts[1], parts[2]);
        }
        s.desc = String::new();
        result.push(s);
    }

    Ok(result)
}

/// CARD database
fn get_card(workdir: &Path, force: bool) -> Result<Vec<DbSeq>> {
    let tarball = workdir.join("card.tar.bz2");
    download("https://card.mcmaster.ca/latest/data", &tarball, force)?;

    eprintln!("  Extracting card.json from tarball...");
    let json_bytes = extract_file_from_tar_bz2(&tarball, "card.json")?;

    // Parse JSON (Latin-1 encoding as in Perl)
    let latin1_string: String = json_bytes.iter().map(|&b| b as char).collect();
    let card: serde_json::Value = serde_json::from_str(&latin1_string)
        .context("Failed to parse card.json")?;

    let mut seqs = Vec::new();

    if let Some(obj) = card.as_object() {
        for (_key, gene) in obj {
            if !gene.is_object() {
                continue;
            }

            let model_type = gene.get("model_type").and_then(|v| v.as_str()).unwrap_or("");
            if model_type != "protein homolog model" {
                continue;
            }

            let id = gene.get("model_name").and_then(|v| v.as_str()).unwrap_or("");
            if id.is_empty() {
                continue;
            }

            // Get DNA sequence
            let dna = gene
                .get("model_sequences")
                .and_then(|ms| ms.get("sequence"))
                .and_then(|s| s.as_object())
                .and_then(|s| s.values().next())
                .and_then(|v| v.get("dna_sequence"));

            let dna = match dna {
                Some(d) => d,
                None => continue,
            };

            let sequence = dna.get("sequence").and_then(|v| v.as_str()).unwrap_or("");
            let accession = dna.get("accession").and_then(|v| v.as_str()).unwrap_or("");
            let fmin = dna.get("fmin").and_then(|v| v.as_i64()).unwrap_or(0);
            let fmax = dna.get("fmax").and_then(|v| v.as_i64()).unwrap_or(0);
            let strand = dna.get("strand").and_then(|v| v.as_str()).unwrap_or("+");

            let (start, stop) = if strand == "-" {
                (fmax, fmin)
            } else {
                (fmin, fmax)
            };

            // Get antibiotic classes
            let mut abx = Vec::new();
            if let Some(aro_cat) = gene.get("ARO_category").and_then(|v| v.as_object()) {
                for (_cat_key, cat) in aro_cat {
                    let class_name = cat.get("category_aro_class_name").and_then(|v| v.as_str()).unwrap_or("");
                    if class_name == "Drug Class" {
                        let mut abx_name = cat.get("category_aro_name").and_then(|v| v.as_str()).unwrap_or("").to_string();
                        abx_name = abx_name.replace(" antibiotic", "").replace(' ', "_");
                        abx.push(abx_name);
                    }
                }
            }

            let desc = gene.get("ARO_description").and_then(|v| v.as_str()).unwrap_or(id);
            let id_clean = id.replace(' ', "_");

            seqs.push(DbSeq {
                id: id_clean,
                acc: format!("{}:{}-{}", accession, start, stop),
                desc: desc.to_string(),
                seq: sequence.to_uppercase(),
                abx,
            });
        }
    }

    eprintln!("  card: {} sequences", seqs.len());
    Ok(seqs)
}

/// NCBI AMRFinderPlus database
fn get_ncbi(workdir: &Path, force: bool) -> Result<Vec<DbSeq>> {
    let afp = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest";
    let fna = workdir.join("amr_cds.ffn");
    let tsv = workdir.join("amr_cds.tsv");

    download(&format!("{}/AMR_CDS.fa", afp), &fna, force)?;
    download(&format!("{}/ReferenceGeneCatalog.txt", afp), &tsv, force)?;

    // Load TSV: key by refseq_nucleotide_accession (column 10, 0-indexed)
    let anno = load_tabular(&tsv, 10)?;
    eprintln!("  ncbi: {} reference records", anno.len());

    let content = fs::read_to_string(&fna)?;
    let raw = parse_fasta(&content);

    let mut seqs = Vec::new();
    for rec in raw {
        // Parse ID: "WP_061158039.1|NG_050200.1|1|1|blaTEM-156|blaTEM|class_A_beta-lactamase_TEM-156"
        // Or old format: "1000909371|WP_061158039.1|NG_050200|1|1|blaTEM-156|blaTEM|class_A_beta-lactamase_TEM-156"
        let fields: Vec<&str> = rec.id.split('|').collect();

        let (acc, fp, fn_, gene, prod) = if fields.len() >= 7 {
            // New format: pi|acc|fp|fn|gene|fam|prod
            (fields[1], fields[2], fields[3], fields[4], fields[6])
        } else if fields.len() >= 8 {
            // Old format: gi|pi|acc|fp|fn|gene|fam|prod
            (fields[2], fields[3], fields[4], fields[5], fields[7])
        } else {
            continue;
        };

        // Skip fusion genes
        let fp: i32 = fp.parse().unwrap_or(0);
        let fn_: i32 = fn_.parse().unwrap_or(0);
        if fp != 1 || fn_ != 1 {
            continue;
        }

        // Ensure accession has version
        let mut acc_v = acc.to_string();
        if !acc_v.contains('.') {
            acc_v.push_str(".1");
        }

        // Look up in TSV
        let t = match anno.get(&acc_v) {
            Some(t) => t,
            None => continue,
        };

        // Only keep core AMR genes
        let scope = t.get("scope").map(|s| s.as_str()).unwrap_or("");
        let gtype = t.get("type").map(|s| s.as_str()).unwrap_or("");
        let subtype = t.get("subtype").map(|s| s.as_str()).unwrap_or("");
        if scope != "core" || gtype != "AMR" || subtype != "AMR" {
            continue;
        }

        let gene_name = gene.to_string();
        let mut product = prod.replace('_', " ");
        if product.is_empty() {
            product = gene_name.clone();
        }

        let refseq_acc = t.get("refseq_nucleotide_accession").cloned().unwrap_or_default();
        let subclass = t.get("subclass").cloned().unwrap_or_default();
        let abx: Vec<String> = if subclass.is_empty() {
            Vec::new()
        } else {
            subclass.split('/').map(|s| s.to_string()).collect()
        };

        let s = rec.seq.to_uppercase();
        let sanitized: String = s.chars().map(|c| if "AGTCN".contains(c) { c } else { 'N' }).collect();

        seqs.push(DbSeq {
            id: gene_name,
            acc: refseq_acc,
            desc: product,
            seq: sanitized,
            abx,
        });
    }

    eprintln!("  ncbi: {} sequences", seqs.len());
    Ok(seqs)
}

/// VFDB database
fn get_vfdb(workdir: &Path, force: bool) -> Result<Vec<DbSeq>> {
    let gz = workdir.join("vfdb.fa.gz");
    download("http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz", &gz, force)?;

    let content = decompress_gz(&gz)?;
    fs::write(workdir.join("vfdb.fa"), &content)?;
    let seqs = load_fasta_file(&workdir.join("vfdb.fa"))?;

    let mut result = Vec::new();
    for mut s in seqs {
        // Parse: >VFG000676(gb|AAD32411) (lef) anthrax toxin lethal factor precursor [...]
        // ID = "VFG000676(gb|AAD32411)"
        // DESC = "(lef) anthrax toxin lethal factor precursor [Anthrax toxin (VF0142)] [Bacillus anthracis str. Sterne]"

        // Extract accession from ID
        if let Some(start) = s.id.find('|') {
            if let Some(end) = s.id[start + 1..].find(')') {
                let acc = &s.id[start + 1..start + 1 + end];
                s.acc = acc.to_string();
            }
        }

        // Extract gene from description: "(lef) anthrax toxin..."
        if let Some(close) = s.desc.find(')') {
            if s.desc.starts_with('(') {
                let gene = &s.desc[1..close];
                if !gene.is_empty() {
                    s.id = gene.to_string();
                }
            }
        }

        result.push(s);
    }

    Ok(result)
}

/// E. coli virulence factors (ecoli_vf)
fn get_ecoli_vf(workdir: &Path, force: bool) -> Result<Vec<DbSeq>> {
    let fasta = workdir.join("ecoli_vf.ffn");
    download(
        "https://github.com/phac-nml/ecoli_vf/raw/master/data/repaired_ecoli_vfs_shortnames.ffn",
        &fasta,
        force,
    )?;

    let seqs = load_fasta_file(&fasta)?;

    let mut result = Vec::new();
    for mut s in seqs {
        // Parse: >VFG000748(gi:2865308) (espF) EspF [EspF (VF0182)] [Escherichia coli O127:H6 str. E2348/69]
        // ID = "VFG000748(gi:2865308)"
        // DESC = "(espF) EspF [EspF (VF0182)] [Escherichia coli O127:H6 str. E2348/69]"

        // Extract accession from ID
        let id_part = s.id.clone();
        if let Some(start) = id_part.find('(') {
            if let Some(end) = id_part.find(')') {
                let inner = &id_part[start + 1..end];
                // inner could be "gi:2865308"
                if let Some(colon) = inner.find(':') {
                    s.acc = inner[colon + 1..].to_string();
                } else {
                    s.acc = inner.to_string();
                }
            }
        }

        // Extract gene from desc and remove strain info
        let desc = s.desc.clone();
        // Remove strain name at end: [Escherichia coli...]
        let desc = regex_remove_strain_brackets(&desc);

        // Extract gene: "(espF) EspF" -> gene = "espF", desc = "EspF"
        if desc.starts_with('(') {
            if let Some(close) = desc.find(')') {
                let gene = &desc[1..close];
                if !gene.is_empty() {
                    s.id = gene.to_string();
                }
                let rest = desc[close + 1..].trim().to_string();
                s.desc = rest;
            }
        } else {
            s.desc = desc;
        }

        result.push(s);
    }

    Ok(result)
}

/// Remove strain name brackets from description: "text [strain info]" -> "text"
fn regex_remove_strain_brackets(desc: &str) -> String {
    // Remove all [...] at the end
    let mut result = desc.to_string();
    while result.ends_with(']') {
        if let Some(open) = result.rfind('[') {
            result = result[..open].trim_end().to_string();
        } else {
            break;
        }
    }
    result
}

/// UPEC/ExPEC virulence factors
fn get_upec_expec_vf(workdir: &Path, force: bool) -> Result<Vec<DbSeq>> {
    let tsv_path = workdir.join("upec_expec_vf.tsv");
    download(
        "https://raw.githubusercontent.com/FordeGenomics/ST167_Code/refs/heads/main/UPEC-ExPEC_VF/UPEC_ExPEC_VF.tsv",
        &tsv_path,
        force,
    )?;

    let content = fs::read_to_string(&tsv_path)?;
    let mut headers: Vec<String> = Vec::new();
    let mut seqs = Vec::new();

    for line in content.lines() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();

        if headers.is_empty() {
            headers = cols.iter().map(|s| s.to_string()).collect();
            continue;
        }

        let row: HashMap<String, String> = headers
            .iter()
            .zip(cols.iter())
            .map(|(h, c)| (h.clone(), c.to_string()))
            .collect();

        let gene = row.get("Gene name").cloned().unwrap_or_default();
        if gene.is_empty() {
            continue;
        }

        let acc = row.get("Accession/Source").cloned().unwrap_or_default();
        let begin = row.get("Begin").cloned().unwrap_or_default();
        let end = row.get("End").cloned().unwrap_or_default();
        let desc = row.get("Description").cloned().unwrap_or_default();
        let class = row.get("Class").cloned().unwrap_or_default();
        let sequence = row.get("Sequence").cloned().unwrap_or_default();

        if sequence.is_empty() {
            continue;
        }

        seqs.push(DbSeq {
            id: gene,
            acc: format!("{}:{}-{}", acc, begin, end),
            desc,
            seq: sequence.to_uppercase(),
            abx: if class.is_empty() { vec![] } else { vec![class] },
        });
    }

    eprintln!("  upec_expec_vf: {} sequences", seqs.len());
    Ok(seqs)
}

/// EcOH serotyping database
fn get_ecoh(workdir: &Path, force: bool) -> Result<Vec<DbSeq>> {
    let fasta = workdir.join("EcOH.fa");
    download(
        "https://raw.githubusercontent.com/katholt/srst2/master/data/EcOH.fasta",
        &fasta,
        force,
    )?;

    let seqs = load_fasta_file(&fasta)?;

    let mut result = Vec::new();
    for mut s in seqs {
        // Parse: >1__fliC__fliC-H1__1 AB028471.1;flagellin;H1
        // ID = "1__fliC__fliC-H1__1"
        // DESC = "AB028471.1;flagellin;H1"

        let id_str = s.id.clone();
        let id_parts: Vec<&str> = id_str.split("__").collect();
        if id_parts.len() >= 3 {
            s.id = id_parts[2].to_string();
        }

        // Parse desc: "AB028471.1;flagellin;H1"
        let desc_parts: Vec<&str> = s.desc.split(';').collect();
        if !desc_parts.is_empty() {
            s.acc = desc_parts[0].to_string();
            s.desc = desc_parts[1..].join(" ");
        }

        result.push(s);
    }

    Ok(result)
}

/// BacMet2 database
fn get_bacmet2(workdir: &Path, force: bool) -> Result<Vec<DbSeq>> {
    let fasta = workdir.join("bacmet2.fa");
    download(
        "http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta",
        &fasta,
        force,
    )?;

    let seqs = load_fasta_file(&fasta)?;

    let mut result = Vec::new();
    for mut s in seqs {
        // Parse: >BAC0098|ctpC|sp|P0A502|CTPC_MYCTU Probable manganese/zinc-exporting
        let id_str = s.id.clone();
        let parts: Vec<&str> = id_str.split('|').collect();
        if parts.len() >= 4 {
            s.id = format!("{}-{}", parts[1], parts[0]);
            s.acc = format!("{}:{}", parts[2], parts[3]);
        }

        result.push(s);
    }

    Ok(result)
}

/// Victors database
fn get_victors(workdir: &Path, force: bool) -> Result<Vec<DbSeq>> {
    let ffn = workdir.join("victors.ffn");
    let faa = workdir.join("victors.faa");

    download("http://phidias.us/victors/downloads/gen_downloads.php", &ffn, force)?;
    download("http://phidias.us/victors/downloads/gen_downloads_protein.php", &faa, force)?;

    // Parse FAA file for GI -> (acc, desc) mapping
    let faa_content = fs::read_to_string(&faa)?;
    let mut gi_map: HashMap<String, (String, String)> = HashMap::new();

    for line in faa_content.lines() {
        if !line.starts_with('>') {
            continue;
        }
        // >gi|115534244|ref|YP_783826.1| hypothetical protein pCJ01p4 [Campylobacter jejuni]
        let header = &line[1..];
        if let Some(rest) = header.strip_prefix("gi|") {
            let parts: Vec<&str> = rest.splitn(2, '|').collect();
            if parts.len() < 2 {
                continue;
            }
            let gi = parts[0];
            let rest2 = parts[1];
            // rest2 = "ref|YP_783826.1| hypothetical protein ..."
            if let Some(pipe) = rest2.find('|') {
                let _db = rest2[..pipe].to_string();
                // Extract accession after "ref|"
                let after_db = &rest2[pipe + 1..];
                let acc = if let Some(pipe2) = after_db.find('|') {
                    after_db[..pipe2].to_string()
                } else {
                    after_db.split_whitespace().next().unwrap_or("").to_string()
                };
                let desc = after_db
                    .split('|')
                    .nth(1)
                    .unwrap_or("")
                    .trim()
                    .to_string();
                gi_map.insert(gi.to_string(), (acc, desc));
            }
        }
    }

    eprintln!("  victors: {} protein records", gi_map.len());

    // Parse FFN file
    let seqs = load_fasta_file(&ffn)?;

    let mut result = Vec::new();
    for mut s in seqs {
        // Parse: >gi|115534241:2616-3152 Campylobacter jejuni plasmid pCJ01, complete sequence
        // Extract GI from ID
        if let Some(rest) = s.id.strip_prefix("gi|") {
            let parts: Vec<&str> = rest.splitn(2, ':').collect();
            if !parts.is_empty() {
                let gi = parts[0];
                let coords = parts.get(1).unwrap_or(&"");
                if let Some((acc, desc)) = gi_map.get(gi) {
                    s.acc = acc.clone();
                    s.desc = desc.clone();
                } else {
                    s.acc = format!("gi|{}:{}", gi, coords);
                    s.desc = "hypothetical protein".to_string();
                }
            }
        }

        result.push(s);
    }

    Ok(result)
}

// ===========================================================================
// Main entry point
// ===========================================================================

/// Download and set up a database
pub fn download_database(
    db_name: &str,
    datadir: &Path,
    force: bool,
    blast_paths: &BlastPaths,
) -> Result<()> {
    // Validate database name
    if !AVAILABLE_DATABASES.contains(&db_name) {
        anyhow::bail!(
            "Unknown database '{}'. Available: {}",
            db_name,
            AVAILABLE_DATABASES.join(", ")
        );
    }

    // Create database directory
    let db_dir = datadir.join(db_name);
    if db_dir.exists() && force {
        eprintln!("Removing existing database: {}", db_dir.display());
        fs::remove_dir_all(&db_dir)?;
    }
    fs::create_dir_all(&db_dir)?;

    // Create temp/src directory
    let src_dir = db_dir.join("src");
    fs::create_dir_all(&src_dir)?;

    eprintln!("Setting up '{}' in '{}'", db_name, db_dir.display());

    // Download and parse the database
    let seqs = match db_name {
        "resfinder" => get_resfinder(&src_dir, force)?,
        "plasmidfinder" => get_plasmidfinder(&src_dir, force)?,
        "megares" => get_megares(&src_dir, force)?,
        "argannot" => get_argannot(&src_dir, force)?,
        "card" => get_card(&src_dir, force)?,
        "ncbi" => get_ncbi(&src_dir, force)?,
        "vfdb" => get_vfdb(&src_dir, force)?,
        "ecoli_vf" => get_ecoli_vf(&src_dir, force)?,
        "upec_expec_vf" => get_upec_expec_vf(&src_dir, force)?,
        "ecoh" => get_ecoh(&src_dir, force)?,
        "bacmet2" => get_bacmet2(&src_dir, force)?,
        "victors" => get_victors(&src_dir, force)?,
        _ => anyhow::bail!("Unknown database: {}", db_name),
    };

    if seqs.is_empty() {
        anyhow::bail!("No sequences found for database '{}'", db_name);
    }

    // Deduplicate
    let seqs = dedupe_seq(seqs);

    // Sort by ID
    let mut seqs = seqs;
    seqs.sort_by(|a, b| a.id.cmp(&b.id));

    // Save sequences file
    let sequences_path = db_dir.join("sequences");
    save_fasta(&sequences_path, db_name, &seqs)?;

    // Build BLAST index
    eprintln!("Formatting BLASTN database: {}", sequences_path.display());
    let db_type = if seqs.iter().any(|s| !is_dna_seq(&s.seq)) {
        "prot"
    } else {
        "nucl"
    };
    eprintln!("  Database type: {}", db_type);
    blast::make_blast_db(blast_paths, &sequences_path, db_type)?;

    eprintln!("Done. {} sequences indexed.", seqs.len());
    Ok(())
}
