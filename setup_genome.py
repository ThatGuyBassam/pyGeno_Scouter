"""
setup_genome.py
===============
One-time genome import for pyGeno Scouter.
Run this script ONCE inside the pyGeno conda environment:

    conda activate pygeno_env
    python C:\pyGeno_Scouter\setup_genome.py

This script creates a pyGeno genome package pointing to Ensembl FTP URLs.
pyGeno will download the chromosome files automatically during import.

Total download: ~3GB for full genome. Time: 1-3 hours depending on connection.
You only need to run this once.
"""

import os
import sys

# ─── PATHS ───────────────────────────────────────────────────────────────────

PACKAGES_DIR = r"C:\pyGeno_Scouter\genome_packages"
os.makedirs(PACKAGES_DIR, exist_ok=True)

# ─── ENSEMBL RELEASE 78 URLS ─────────────────────────────────────────────────

ENSEMBL_FASTA_BASE = "ftp://ftp.ensembl.org/pub/release-78/fasta/homo_sapiens/dna"
ENSEMBL_GTF        = "ftp://ftp.ensembl.org/pub/release-78/gtf/homo_sapiens/Homo_sapiens.GRCh38.78.gtf.gz"

# All human chromosomes in GRCh38 release 78
CHROMOSOMES = {
    "1":  f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz",
    "2":  f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz",
    "3":  f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz",
    "4":  f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz",
    "5":  f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz",
    "6":  f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz",
    "7":  f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz",
    "8":  f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.8.fa.gz",
    "9":  f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.9.fa.gz",
    "10": f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz",
    "11": f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz",
    "12": f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz",
    "13": f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz",
    "14": f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz",
    "15": f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz",
    "16": f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.16.fa.gz",
    "17": f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz",
    "18": f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.18.fa.gz",
    "19": f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz",
    "20": f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz",
    "21": f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz",
    "22": f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz",
    "X":  f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz",
    "Y":  f"{ENSEMBL_FASTA_BASE}/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz",
}

# ─── STEP 1: Create Y-only test package ──────────────────────────────────────

def create_y_only_package():
    pkg_dir = os.path.join(PACKAGES_DIR, "GRCh38.78_Y-only")
    os.makedirs(pkg_dir, exist_ok=True)

    manifest = f"""[package_infos]
description = GRCh38.78 Y chromosome only — test package
maintainer = pyGeno Scouter setup
maintainer_contact = N/A
version = GRCh38.78_Y-only

[genome]
species = Homo_sapiens
name = GRCh38.78_Y-only
source = http://www.ensembl.org

[chromosome_files]
Y = {CHROMOSOMES['Y']}

[gene_set]
gtf = {ENSEMBL_GTF}
"""
    with open(os.path.join(pkg_dir, "manifest.ini"), "w") as f:
        f.write(manifest)

    return pkg_dir


def create_full_package():
    pkg_dir = os.path.join(PACKAGES_DIR, "GRCh38.78")
    os.makedirs(pkg_dir, exist_ok=True)

    chrom_lines = "\n".join(
        f"{chrom} = {url}" for chrom, url in CHROMOSOMES.items()
    )

    manifest = f"""[package_infos]
description = GRCh38.78 full human genome
maintainer = pyGeno Scouter setup
maintainer_contact = N/A
version = GRCh38.78

[genome]
species = Homo_sapiens
name = GRCh38.78
source = http://www.ensembl.org

[chromosome_files]
{chrom_lines}

[gene_set]
gtf = {ENSEMBL_GTF}
"""
    with open(os.path.join(pkg_dir, "manifest.ini"), "w") as f:
        f.write(manifest)

    return pkg_dir


# ─── MAIN ────────────────────────────────────────────────────────────────────

from pyGeno.importation.Genomes import importGenome

print("=" * 60)
print("pyGeno Genome Setup")
print("=" * 60)

# Step 1 — Y-only test
print("\n[Step 1] Creating Y-chromosome test package...")
y_pkg = create_y_only_package()
print(f"         Package folder: {y_pkg}")

print("\n[Step 1] Importing Y-only genome (downloads from Ensembl FTP)...")
print("         This may take a few minutes.\n")
try:
    importGenome(y_pkg)
    print("[Step 1] Y-only genome imported successfully.")
except Exception as e:
    print(f"[Step 1] ERROR: {e}")
    sys.exit(1)

# Step 2 — Verify Y-only
print("\n[Step 2] Verifying Y-only genome...")
try:
    from pyGeno.Genome import Genome
    g = Genome(name="GRCh38.78_Y-only")
    print(f"[Step 2] OK — genome loaded.")
except Exception as e:
    print(f"[Step 2] ERROR: {e}")
    sys.exit(1)

# Step 3 — Full genome
print("\n[Step 3] Creating full GRCh38.78 genome package...")
full_pkg = create_full_package()
print(f"         Package folder: {full_pkg}")

print("\n[Step 3] Importing full genome (downloads ~3GB from Ensembl FTP)...")
print("         This will take 1-3 hours. Do NOT close this terminal.\n")
try:
    importGenome(full_pkg)
    print("\n[Step 3] Full genome imported successfully.")
except Exception as e:
    print(f"\n[Step 3] ERROR: {e}")
    print("If this was a network timeout, run the script again — it may resume.")
    sys.exit(1)

# Step 4 — Final check
print("\n[Step 4] Final verification...")
try:
    from pyGeno.Genome import Genome
    g   = Genome(name="GRCh38.78")
    hbb = g.get("Gene", name="HBB")[0]
    print(f"[Step 4] HBB found on chromosome {hbb.chromosome.number}")
    print(f"         Position: {hbb.start:,} – {hbb.end:,}")
    print("\n" + "=" * 60)
    print("Setup complete. Launch pyGeno Scouter:")
    print("  streamlit run C:\\pyGeno_Scouter\\app.py")
    print("=" * 60)
except Exception as e:
    print(f"[Step 4] ERROR: {e}")
    sys.exit(1)
