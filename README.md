# 🧬 pyGeno Scouter
**Clinical phenotype → genomic coordinates**

A Streamlit interface that maps clinical symptoms and disease names to candidate genes, genomic coordinates, transcript isoforms, and known pathogenic variants. Built on top of [pyGeno](https://github.com/tariqdaouda/pyGeno) by Tariq Daouda (IRIC Montréal), with a curated database oriented around the North African patient population.

---

## Offline / Online split

| Feature | Connection |
|---|---|
| Phenotype Search (curated, HPO, Orphanet) | ✅ Fully offline |
| pyGeno genomic coordinates & sequences | ✅ Fully offline |
| Gene Batch — ClinVar pathogenic variants | 🌐 Requires internet |

The core use case — typing a symptom and getting candidate genes and genomic data — works with no internet connection after initial setup. ClinVar queries in the Gene Batch tab hit the NCBI E-utilities API live.

---

## What it does

You type a clinical phenotype — a symptom, a disease name, a gene name, or a syndrome — and the tool returns:

- Candidate genes associated with that phenotype
- Genomic coordinates (chromosome, position, strand) from Ensembl GRCh38
- Transcript isoforms with exon-by-exon breakdowns
- Protein sequences for each coding isoform
- Known pathogenic variants with clinical descriptions

A second tab (**Gene Batch + ClinVar**) accepts a list of gene names — one per line or comma-separated — and returns pyGeno coordinates plus live ClinVar pathogenic/likely pathogenic variants for each gene, with a one-click CSV export of all variants.

### Search layers

| Layer | Source | Language | Coverage |
|---|---|---|---|
| Curated | Hand-curated, North Africa-focused | French / English | 13 conditions |
| HPO | Human Phenotype Ontology (JAX) | English | 42,553 phenotype terms |
| Orphanet | Orphanet rare disease database | French / English | 4,128 diseases |
| ClinVar | NCBI ClinVar (live) | English | Pathogenic variants by gene |

---

## The problem it solves

pyGeno is a powerful tool for querying personalized genome data programmatically — but it requires Python literacy to use. A clinician or medical student who wants to know the genomic coordinates of HBB, or which genes are associated with recurrent fever in a North African patient, has to write Python code to get that information out of pyGeno. Most do not.

At the same time, the major disease-to-gene databases — OMIM, HPO, Orphanet — exist separately from each other and from pyGeno. A clinician looking up a patient with hemolytic anemia would need to consult OMIM for gene-disease links, HPO for phenotype mappings, Orphanet for rare disease coverage, and then separately query a genome browser for coordinates. None of these databases talk to each other natively.

pyGeno Scouter connects all layers in a single interface:

```
Clinical language  →  Candidate genes  →  Live genomic data  →  Pathogenic variants
(HPO / Orphanet)      (OMIM / curated)     (pyGeno / Ensembl)     (ClinVar / NCBI)
```

It also addresses a specific gap in existing tools: the North African patient population is underrepresented in most genomic databases. The curated database in this project is built around conditions with elevated prevalence in Moroccan and North African patients — thalassemia, G6PD deficiency, FMF, CFTR mutations — with clinical notes and variant data reflecting that specific population.

---

## Who is it for

**Medical students and residents** who want to connect a clinical presentation to its genetic basis without writing code.

**Researchers using pyGeno** who want a fast interface for exploratory queries — checking gene positions, isoform structures, protein sequences, or known pathogenic variants — without writing a script for every lookup.

**Clinicians in North African and Mediterranean settings** where conditions like beta-thalassemia, sickle cell disease, G6PD deficiency, and FMF are clinically common but often underserved by tools calibrated to European reference populations.

**Anyone building on pyGeno** who wants to see how to bridge the Python 3.6 / Python 3.8+ compatibility gap using a subprocess architecture.

---

## Requirements

- Python 3.11 (for Streamlit)
- Miniconda (for the Python 3.6 pyGeno environment)
- ~4GB disk space for the genome

---

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/Thatguybassam/pygeno-scouter
cd pygeno-scouter
```

### 2. Install Streamlit dependencies
```bash
pip install streamlit pandas
```

### 3. Download the data files

The three database files are not included in the repo. Download them and place in `data/`:

| File | Source |
|---|---|
| `hp.obo` | https://hpo.jax.org/app/data/ontology |
| `phenotype_to_genes.txt` | https://hpo.jax.org/app/data/annotations |
| `en_product6.xml` | https://www.orphadata.com/genes/ |

```
data/
├── hp.obo
├── phenotype_to_genes.txt
└── en_product6.xml
```

### 4. Set up the pyGeno environment

pyGeno requires Python 3.6. Use Miniconda:

```bash
conda create -n pygeno_env python=3.6
conda activate pygeno_env
pip install pyGeno rabaDB
```

Patch the version check — pyGeno's version guard incorrectly rejects Python 3. Open `<miniconda_path>/envs/pygeno_env/Lib/site-packages/pyGeno/configuration.py` and comment out:

```python
# if not checkPythonVersion():
#     raise PythonVersionError(...)
```

Import the genome — one-time, downloads ~3GB from Ensembl FTP, takes 1–3 hours:

```bash
conda activate pygeno_env
python setup_genome.py
```

### 5. Configure the pyGeno path

Open `app.py` and update:

```python
PYGENO_PYTHON = r"C:\Users\<your_username>\miniconda3\envs\pygeno_env\python.exe"
```

---

## Running the app

```bash
streamlit run app.py
```

Or double-click `launch_scouter.bat` on Windows. Opens at http://localhost:8501.

---

## Project structure

```
pygeno-scouter/
├── app.py               # Streamlit interface
├── clinvar.py           # ClinVar API module (NCBI E-utilities)
├── phenotype_db.py      # Curated North Africa disease database
├── data_loader.py       # HPO and Orphanet parsers + search functions
├── pygeno_query.py      # pyGeno backend (subprocess, Python 3.6)
├── setup_genome.py      # One-time genome import script
├── launch_scouter.bat   # Windows launcher
├── SOURCES.txt          # Full data source citations
├── .gitignore
└── data/                # Local only — not committed
    ├── hp.obo
    ├── phenotype_to_genes.txt
    └── en_product6.xml
```

---

## Architecture

The app runs two Python environments simultaneously. Streamlit requires Python 3.8+. pyGeno requires Python 3.6. They cannot share an environment.

`pygeno_query.py` runs as a subprocess inside the Python 3.6 conda environment. The Streamlit app calls it like a command-line tool, passes a gene name as an argument, and reads back a JSON response. This keeps the environments fully isolated while allowing live genomic queries.

ClinVar queries use only Python stdlib (`urllib`, `json`) — no extra dependencies required. Results are cached for one hour per session to avoid redundant API calls.

The search pipeline:

```
User query
    ↓
phenotype_db.py     → curated North Africa results (French/English, offline)
data_loader.py      → HPO + Orphanet results (offline)
    ↓
For each gene found:
pygeno_query.py     → live genomic data via pyGeno (offline)
    ↓
Rendered in Streamlit

Gene Batch tab:
gene list → pygeno_query.py (offline) + clinvar.py (internet) → CSV export
```

---

## Data sources

| Source | Use |
|---|---|
| pyGeno — Daouda et al., F1000Research 2016 | Genomic coordinates, sequences, isoforms |
| Ensembl GRCh38.78 | Reference genome |
| HPO — Köhler et al., NAR 2021 | Phenotype → gene associations |
| Orphanet — CC-BY-4.0 | Rare disease → gene associations |
| ClinVar — NCBI | Pathogenic variant classifications (live) |

Full citations with specific entries in `SOURCES.txt`.

---

## Curated conditions

13 conditions with elevated prevalence or clinical relevance in North African populations:

Anémie Hémolytique · Drépanocytose · Bêta-Thalassémie · Déficience en G6PD · Dystrophies Musculaires (DMD, DYSF, CAPN3) · Neuropathie Périphérique (CMT) · Épilepsie · Fièvre Récurrente / FMF · Déficit Immunitaire · Phénylcétonurie · Diabète Monogénique · Mucoviscidose · Cardiomyopathie Hypertrophique · Syndrome du QT Long

---

## Demo

<!-- Add screenshot or GIF here showing Phenotype Search and Gene Batch tabs -->
<!-- Recommended: record a short GIF using ScreenToGif (Windows) or Peek (Linux) -->
<!-- showing: typing "thalassemia" in Phenotype Search, then pasting HBB/G6PD/MEFV in Gene Batch -->

---

## Limitations

- pyGeno is a genome browser, not a diagnostic tool. It contains no clinical decision logic.
- The curated database is a manually maintained lookup table, not a differential diagnosis engine.
- HPO is indexed in English only. French queries work better against Orphanet or the curated layer.
- ClinVar results depend on internet availability. The rest of the app is unaffected if offline.
- Not validated for clinical use. Do not use to make clinical decisions.

---

## License

MIT

pyGeno is the work of Tariq Daouda at IRIC Montréal. Orphanet data: CC-BY-4.0. HPO data: see hpo.jax.org for license terms.
