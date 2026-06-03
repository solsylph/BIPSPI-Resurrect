#!/usr/bin/env bash
# Phase D: download UniRef90 and compile a BLAST protein database.
#
# Runtime: ~5-10 hours total, dominated by the ~30 GB download.
# Disk: ~30 GB compressed -> ~70 GB after gunzip + makeblastdb. You have
# ~17 TB free on /home so this fits comfortably.
#
# Run inside tmux + srun --pty so an SSH disconnect doesn't kill it:
#   tmux new -s uniref90
#   srun --partition=cpu --cpus-per-task=4 --mem=16G --time=12:00:00 --pty bash -l
#   conda activate protein
#   module load BLAST+/2.14.1-gompi-2023a
#   bash tools/download_uniref90_db.sh
#   # Ctrl-b d to detach.  Reattach later with: tmux attach -t uniref90
#
# When done, the BLAST DB lives at:
#   ~/databases/uniref90/uniref90.fasta  (the .pal/.phr/.pin/.psq files share this prefix)
# Use that as psiBlastDB_path in dependencies.cfg.

set -euo pipefail

DB_DIR="${DB_DIR:-$HOME/databases/uniref90}"
mkdir -p "$DB_DIR"
cd "$DB_DIR"

# Verify makeblastdb is on PATH (must be loaded with `module load BLAST+/...`)
if ! command -v makeblastdb >/dev/null 2>&1; then
  echo "ERROR: makeblastdb not on PATH."
  echo "Run: module load BLAST+/2.14.1-gompi-2023a"
  exit 1
fi

# ---------------------------------------------------------------------------
# 1. Download uniref90.fasta.gz from UniProt FTP
# ---------------------------------------------------------------------------
FASTA_GZ=uniref90.fasta.gz
FASTA=uniref90.fasta

if [ ! -f "$FASTA" ] && [ ! -f "$FASTA_GZ" ]; then
  echo "=== Downloading UniRef90 (~30 GB compressed) ==="
  # UniProt's official source. -c resumes on interrupt.
  wget -c https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
fi

# ---------------------------------------------------------------------------
# 2. Decompress
# ---------------------------------------------------------------------------
if [ ! -f "$FASTA" ]; then
  echo "=== Decompressing (~30 min) ==="
  gunzip -v "$FASTA_GZ"
fi

# ---------------------------------------------------------------------------
# 3. Compile BLAST DB
# ---------------------------------------------------------------------------
if [ ! -f "${FASTA}.pal" ] && [ ! -f "${FASTA}.phr" ]; then
  echo "=== Compiling BLAST DB (makeblastdb, ~60-90 min) ==="
  # Minimal flags only: -in / -dbtype / -out.
  # Dropped -parse_seqids (memory) and -hash_index (the two flags interact
  # poorly without -parse_seqids -> "Duplicate seq_ids" error on UniRef90).
  # psiblast queries work fine against the basic pin/phr/psq indexes.
  makeblastdb \
    -in "$FASTA" \
    -dbtype prot \
    -out "$FASTA"
fi

echo ""
echo "=== DONE ==="
echo "BLAST DB compiled at: $DB_DIR/$FASTA"
echo ""
echo "Now edit configFiles/cmdTool/dependencies.cfg to set:"
echo "    psiBlastDB_path   $DB_DIR/$FASTA"
echo ""
echo "Verify it works:"
echo "    echo -e '>test\\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF' > /tmp/test.fa"
echo "    psiblast -query /tmp/test.fa -db $DB_DIR/$FASTA -num_iterations 1 -num_threads 4 -out /tmp/test.psiblast"
echo "    head -50 /tmp/test.psiblast"
