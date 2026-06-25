# abricate_rust

Rust implementation of [ABRicate](https://github.com/tseemann/abricate) — mass screening of contigs for antimicrobial resistance (AMR) and virulence genes. Cross-platform (Windows, Linux, macOS) with integrated [any2fasta](https://github.com/tseemann/any2fasta) and configurable BLAST+ path.

## Features

- **Cross-platform**: Runs natively on Windows, Linux, and macOS (no Perl required)
- **Integrated any2fasta**: Built-in sequence format conversion (FASTA, GenBank, EMBL, GFF, FASTQ, GFA, Clustal, Stockholm, PDB) — no external dependency
- **Configurable BLAST+**: Point to any BLAST+ installation via `--blastdir` or `ABRICATE_BLAST_DIR` environment variable
- **Multiple output formats**: TSV (default), CSV, BED, GFF3, JSON
- **Summary mode**: Generate gene presence/absence matrix across multiple samples
- **Database compatible**: Uses the same `~~~`-separated FASTA header format as the original ABRicate

## Quick Start

### 1. Install BLAST+

Download BLAST+ from [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and extract to a directory.

### 2. Build abricate_rust

```bash
cargo build --release
```

The binary will be at `target/release/abricate` (or `abricate.exe` on Windows).

### 3. Check dependencies

```bash
abricate --blastdir /path/to/ncbi-blast/bin --check
```

### 4. Set up databases

```bash
abricate --blastdir /path/to/ncbi-blast/bin --setupdb
```

### 5. Run detection

```bash
abricate --blastdir /path/to/ncbi-blast/bin --db ncbi genome.fasta
```

## Usage

```
abricate [OPTIONS] <FASTA_FILES>...

Options:
      --list                 List all available databases
      --setupdb              Build BLAST indices for all databases
      --check                Check that all required dependencies are available
      --summary              Generate a summary matrix from abricate reports
      --identity             Use %IDENTITY instead of %COVERAGE in summary mode
      --db <DB>              Database to use [default: ncbi]
      --datadir <DATADIR>    Data directory containing databases [default: auto-detect]
      --blastdir <BLASTDIR>  Directory containing BLAST+ binaries
                             Can also be set via ABRICATE_BLAST_DIR env var
      --minid <MINID>        Minimum DNA %identity [default: 80]
      --mincov <MINCOV>      Minimum DNA %coverage [default: 80]
  -t, --threads <THREADS>    Number of BLAST threads [default: 1]
      --outfmt <OUTFMT>      Output format: tsv, csv, bed, gff, json [default: tsv]
      --noheader             Suppress header row in output
      --nopath               Strip directory paths from FILE column
      --fofn <FOFN>          File containing a list of input files (one per line)
  -h, --help                 Print help
  -V, --version              Print version
```

## BLAST+ Configuration

abricate_rust needs BLAST+ binaries (`blastn`, `blastx`, `makeblastdb`, `blastdbcmd`). You can specify their location in three ways:

1. **Command-line flag** (highest priority):
   ```bash
   abricate --blastdir /path/to/ncbi-blast/bin input.fasta
   ```

2. **Environment variable**:
   ```bash
   export ABRICATE_BLAST_DIR=/path/to/ncbi-blast/bin
   abricate input.fasta
   ```
   
   On Windows:
   ```cmd
   set ABRICATE_BLAST_DIR=D:\code\ncbi-blast-2.17.0+\bin
   abricate input.fasta
   ```

3. **System PATH** (lowest priority): If BLAST+ is already in your PATH, no configuration is needed.

## Database Format

Each database is a subdirectory under the data directory containing a `sequences` FASTA file. FASTA headers use `~~~` as a field separator:

```
>DATABASE~~~GENE~~~ACCESSION~~~RESISTANCE Product description
ATGCGATCGATCGATCGATCG...
```

Example directory structure:
```
db/
├── ncbi/
│   └── sequences          # FASTA file with AMR gene sequences
├── card/
│   └── sequences
├── vfdb/
│   └── sequences
└── resfinder/
    └── sequences
```

After running `--setupdb`, BLAST index files (`.nhr`, `.nin`, `.nsq`) will be created alongside the `sequences` file.

## Output Format

### TSV (default)

```
#FILE    SEQUENCE  START  END  STRAND  GENE    COVERAGE   %COVERAGE  %IDENTITY  DATABASE  ACCESSION  PRODUCT              RESISTANCE
sample.fa  contig_1  201    700  +       blaTEM-1  1-500/500  100.00     100.00     ncbi      TEST001    beta-lactamase TEM-1  BETA-LACTAM
```

### Summary Mode

```
abricate --summary report1.tsv report2.tsv report3.tsv
```

Generates a gene presence/absence matrix:
```
#FILE       blaTEM-1  aac3-Ia  tetA
report1      100.0     100.0    .
report2      100.0       .      .
report3        .       100.0    .
```

## Environment Variables

| Variable | Description |
|----------|-------------|
| `ABRICATE_BLAST_DIR` | Path to BLAST+ binaries directory |
| `ABRICATE_DATADIR` | Path to database directory |

## Building from Source

```bash
git clone https://github.com/fanvanf/abricate_rust.git
cd abricate_rust
cargo build --release
```

## Differences from Original ABRicate

| Feature | Original (Perl) | Rust Version |
|---------|-----------------|--------------|
| Runtime | Perl + external any2fasta | Single native binary |
| Platform | Linux/macOS (limited Windows) | Windows, Linux, macOS |
| any2fasta | External dependency | Integrated in-process |
| BLAST+ path | Must be in PATH | Configurable via flag/env |
| Pipeline | `bash -c` pipe | In-memory with stdin piping |
| Output | TSV, CSV, BED, GFF | TSV, CSV, BED, GFF, JSON |

## License

GPL-2.0, same as the original ABRicate.

## Acknowledgments

- [Torsten Seemann](https://github.com/tseemann) for the original [ABRicate](https://github.com/tseemann/abricate) and [any2fasta](https://github.com/tseemann/any2fasta) tools
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) team
