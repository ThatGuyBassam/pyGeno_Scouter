"""
app.py — pyGeno Scouter
=======================
Clinical phenotype → genomic coordinates interface.
Curated North Africa-focused database + HPO/Orphanet fallback.

Run:
    streamlit run app.py

Requirements (Python 3.11):
    pip install streamlit pandas
"""

import os
import sys
import json
import subprocess
import streamlit as st
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from phenotype_db import PHENOTYPE_DB, search_phenotypes
from data_loader import (
    load_hpo_ontology, load_hpo_gene_map, load_orphanet,
    search_hpo, search_orphanet
)

# ─── CONFIGURATION ──────────────────────────────────────────────────────────

PYGENO_PYTHON = r"C:\Users\GAMER\miniconda3\envs\pygeno_env\python.exe"
PYGENO_SCRIPT = os.path.join(os.path.dirname(__file__), "pygeno_query.py")
GENOME_BUILD  = "GRCh38.78"

# ─── PAGE CONFIG ────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="pyGeno Scouter",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─── CSS ────────────────────────────────────────────────────────────────────

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap');

/* ── Reset & Base ── */
*, *::before, *::after { box-sizing: border-box; }

html, body, .stApp {
    background-color: #0d0f14 !important;
    color: #c8cfe0 !important;
    font-family: 'IBM Plex Sans', sans-serif !important;
}

.block-container {
    padding: 2.5rem 3rem !important;
    max-width: 1080px !important;
}

/* ── Hide Streamlit chrome ── */
#MainMenu { display: none; }
footer { visibility: hidden; }
.stDeployButton { display: none; }

/* ── Sidebar ── */
[data-testid="stSidebar"] {
    background-color: #0a0c10 !important;
    border-right: 1px solid #1e2230 !important;
}
[data-testid="stSidebar"] * {
    color: #8090a8 !important;
    font-family: 'IBM Plex Sans', sans-serif !important;
    font-size: 13px !important;
}
[data-testid="stSidebar"] h2 {
    color: #c8cfe0 !important;
    font-size: 16px !important;
    font-weight: 600 !important;
    letter-spacing: 0.02em !important;
}
[data-testid="stSidebar"] .stSuccess {
    background: #0f1f10 !important;
    border: 1px solid #1e3a20 !important;
    color: #4a9e5c !important;
    border-radius: 6px !important;
    font-size: 12px !important;
}
[data-testid="stSidebar"] .stWarning {
    background: #1a1200 !important;
    border: 1px solid #3a2800 !important;
    color: #8a7020 !important;
    border-radius: 6px !important;
    font-size: 12px !important;
}

/* ── Title ── */
h1 {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 26px !important;
    font-weight: 600 !important;
    color: #e8eeff !important;
    letter-spacing: -0.01em !important;
    margin-bottom: 4px !important;
}

/* ── Search input ── */
[data-testid="stTextInput"] input {
    background: #131720 !important;
    border: 1px solid #2a3050 !important;
    border-radius: 8px !important;
    color: #c8cfe0 !important;
    font-family: 'IBM Plex Sans', sans-serif !important;
    font-size: 15px !important;
    padding: 12px 16px !important;
    transition: border-color 0.2s !important;
}
[data-testid="stTextInput"] input:focus {
    border-color: #4a6cf7 !important;
    box-shadow: 0 0 0 3px rgba(74,108,247,0.12) !important;
    outline: none !important;
}
[data-testid="stTextInput"] input::placeholder {
    color: #445068 !important;
}

/* ── Primary button ── */
.stButton > button[kind="primary"] {
    background: #4a6cf7 !important;
    border: none !important;
    border-radius: 8px !important;
    color: #ffffff !important;
    font-family: 'IBM Plex Sans', sans-serif !important;
    font-weight: 600 !important;
    font-size: 14px !important;
    padding: 12px 20px !important;
    transition: background 0.2s !important;
}
.stButton > button[kind="primary"]:hover {
    background: #3a5ce7 !important;
}
.stButton > button {
    background: #1a1f30 !important;
    border: 1px solid #2a3050 !important;
    border-radius: 8px !important;
    color: #8090b8 !important;
    font-family: 'IBM Plex Sans', sans-serif !important;
}

/* ── Expander ── */
[data-testid="stExpander"] {
    background: #111420 !important;
    border: 1px solid #1e2438 !important;
    border-radius: 10px !important;
    margin: 6px 0 !important;
    overflow: hidden !important;
}
[data-testid="stExpander"] summary {
    background: #111420 !important;
    color: #c8cfe0 !important;
    font-family: 'IBM Plex Sans', sans-serif !important;
    font-size: 14px !important;
    font-weight: 500 !important;
    padding: 14px 18px !important;
}
[data-testid="stExpander"] summary:hover {
    background: #161b2e !important;
}
[data-testid="stExpander"] > div > div {
    background: #0e1118 !important;
    padding: 16px 18px !important;
}

/* ── Tabs ── */
[data-testid="stTabs"] [data-baseweb="tab-list"] {
    background: transparent !important;
    border-bottom: 1px solid #1e2438 !important;
    gap: 0 !important;
}
[data-testid="stTabs"] [data-baseweb="tab"] {
    background: transparent !important;
    color: #5060a0 !important;
    font-family: 'IBM Plex Sans', sans-serif !important;
    font-size: 13px !important;
    font-weight: 500 !important;
    padding: 8px 16px !important;
    border-bottom: 2px solid transparent !important;
}
[data-testid="stTabs"] [aria-selected="true"] {
    color: #8090e0 !important;
    border-bottom-color: #4a6cf7 !important;
}
[data-testid="stTabs"] [data-baseweb="tab-panel"] {
    background: transparent !important;
    padding: 16px 0 !important;
}

/* ── Metrics ── */
[data-testid="stMetric"] {
    background: #131720 !important;
    border: 1px solid #1e2438 !important;
    border-radius: 8px !important;
    padding: 12px 16px !important;
}
[data-testid="stMetricLabel"] {
    color: #506080 !important;
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 11px !important;
    text-transform: uppercase !important;
    letter-spacing: 0.06em !important;
}
[data-testid="stMetricValue"] {
    color: #c8cfe0 !important;
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 18px !important;
    font-weight: 600 !important;
}

/* ── Dataframe ── */
[data-testid="stDataFrame"] {
    border: 1px solid #1e2438 !important;
    border-radius: 8px !important;
    overflow: hidden !important;
}

/* ── Info / Warning ── */
.stInfo {
    background: #0d1525 !important;
    border: 1px solid #1e3060 !important;
    border-radius: 8px !important;
    color: #6080c0 !important;
}
.stError {
    background: #1a0d0d !important;
    border: 1px solid #401010 !important;
    color: #c06060 !important;
    border-radius: 8px !important;
}

/* ── Custom components ── */
.gene-chip {
    display: inline-block;
    background: #141c38;
    color: #7090e8;
    border: 1px solid #2a3870;
    border-radius: 20px;
    padding: 3px 12px;
    margin: 3px 2px;
    font-size: 12px;
    font-weight: 600;
    font-family: 'IBM Plex Mono', monospace;
    letter-spacing: 0.02em;
}
.snp-row {
    background: #0e1220;
    border-left: 2px solid #4a6cf7;
    border-radius: 0 6px 6px 0;
    padding: 8px 12px;
    margin: 4px 0;
    font-size: 13px;
    color: #a0b0d0;
    font-family: 'IBM Plex Sans', sans-serif;
}
.snp-row code {
    font-family: 'IBM Plex Mono', monospace;
    color: #7090e8;
    font-size: 12px;
}
.clinical-box {
    background: #0e1525;
    border-left: 3px solid #4a6cf7;
    border-radius: 0 8px 8px 0;
    padding: 14px 18px;
    margin: 12px 0;
    color: #8090c0;
    font-size: 13px;
    line-height: 1.7;
}
.genomic-pos {
    display: inline-block;
    font-family: 'IBM Plex Mono', monospace;
    font-size: 12px;
    color: #c08040;
    background: #1a1408;
    border: 1px solid #2a2010;
    border-radius: 6px;
    padding: 5px 12px;
    margin: 8px 0;
}
.section-label {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 10px;
    font-weight: 600;
    letter-spacing: 0.12em;
    color: #3a4868;
    text-transform: uppercase;
    margin: 18px 0 8px 0;
}
.badge {
    display: inline-block;
    border-radius: 4px;
    padding: 2px 8px;
    font-size: 11px;
    font-weight: 600;
    font-family: 'IBM Plex Mono', monospace;
    letter-spacing: 0.04em;
}
.badge-curated  { background: #0f1f10; color: #4a8e5a; border: 1px solid #1e3a20; }
.badge-hpo      { background: #101828; color: #4a6ec8; border: 1px solid #1e3060; }
.badge-orphanet { background: #180f28; color: #8a5ec8; border: 1px solid #301a50; }
.badge-icd      { background: #180f28; color: #8a5ec8; border: 1px solid #301a50; }
.tip-box {
    background: #0e1220;
    border: 1px solid #1e2438;
    border-radius: 8px;
    padding: 14px 18px;
    color: #506080;
    font-size: 13px;
    line-height: 1.6;
}
.divider {
    border: none;
    border-top: 1px solid #1a1e2e;
    margin: 20px 0;
}
</style>
""", unsafe_allow_html=True)

# ─── LOAD DATA ───────────────────────────────────────────────────────────────

@st.cache_resource(show_spinner="Loading genomic databases...")
def load_all_data():
    term_to_id, id_to_name    = load_hpo_ontology()
    hp_to_genes, gene_to_hps  = load_hpo_gene_map()
    orphanet_diseases         = load_orphanet()
    return term_to_id, id_to_name, hp_to_genes, gene_to_hps, orphanet_diseases

term_to_id, id_to_name, hp_to_genes, gene_to_hps, orphanet_diseases = load_all_data()

# ─── PYGENO ──────────────────────────────────────────────────────────────────

@st.cache_data(ttl=3600, show_spinner=False)
def run_pygeno_query(gene_name: str) -> dict:
    if not os.path.exists(PYGENO_PYTHON):
        return {"error": "not_configured"}
    try:
        result = subprocess.run(
            [PYGENO_PYTHON, PYGENO_SCRIPT, gene_name, GENOME_BUILD],
            capture_output=True, text=True, timeout=60
        )
        if result.returncode != 0:
            return {"error": result.stderr.strip()}
        return json.loads(result.stdout)
    except subprocess.TimeoutExpired:
        return {"error": "Query timed out (>60s). Is the genome imported?"}
    except Exception as e:
        return {"error": str(e)}


def pygeno_available():
    return os.path.exists(PYGENO_PYTHON)


def protein_blocks(seq, block=10, per_line=6):
    if not seq:
        return ""
    lines = []
    for ls in range(0, len(seq), block * per_line):
        chunk = seq[ls:ls + block * per_line]
        blocks = [chunk[i:i+block] for i in range(0, len(chunk), block)]
        lines.append(f"{ls+1:>5}  " + "  ".join(blocks))
    return "\n".join(lines)


def show_pygeno_section(genes: list):
    if not pygeno_available():
        st.markdown(
            '<div class="tip-box">⚙ pyGeno not configured — '
            'set <code>PYGENO_PYTHON</code> in app.py and import the genome '
            'to enable live sequence data.</div>',
            unsafe_allow_html=True
        )
        return

    tabs = st.tabs(genes[:6])
    for tab, gene_name in zip(tabs, genes[:6]):
        with tab:
            with st.spinner(f"Querying pyGeno — {gene_name}..."):
                data = run_pygeno_query(gene_name)

            if "error" in data:
                if data["error"] == "not_configured":
                    st.info("Set PYGENO_PYTHON in app.py to enable genomic queries.")
                else:
                    st.error(f"pyGeno: {data['error']}")
                continue

            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Chromosome", data["chromosome"])
            c2.metric("Length", f"{data['length_bp']:,} bp")
            c3.metric("Transcripts", data["n_transcripts"])
            c4.metric("Coding isoforms", data["n_coding_transcripts"])

            strand = "Forward (+)" if data["strand"] == 1 else "Reverse (−)"
            st.markdown(
                f'<div class="genomic-pos">chr{data["chromosome"]}:'
                f'{data["start"]:,}–{data["end"]:,} &nbsp;·&nbsp; {strand}'
                f' &nbsp;·&nbsp; {data["ensembl_id"]}</div>',
                unsafe_allow_html=True
            )

            for iso in data.get("isoforms", []):
                with st.expander(
                    f'{iso["transcript_id"]} — {iso["protein_length"]} aa '
                    f'· {iso["coding_exons"]} coding exons'
                ):
                    if iso["exons"]:
                        df = pd.DataFrame(iso["exons"])
                        df.columns = ["Start", "End", "CDS (bp)", "5' CDS"]
                        df.index = [f"Exon {i+1}" for i in range(len(df))]
                        df["Start"] = df["Start"].apply(lambda x: f"{x:,}")
                        df["End"]   = df["End"].apply(lambda x: f"{x:,}")
                        st.dataframe(df, use_container_width=True)
                    if iso.get("protein_seq_60"):
                        st.code(
                            protein_blocks(iso["protein_seq_60"]) +
                            (f"\n... [{iso['protein_length'] - 60} more amino acids]"
                             if iso["protein_length"] > 60 else ""),
                            language=None
                        )

# ─── RENDERERS ───────────────────────────────────────────────────────────────

def render_curated_result(condition_key, condition):
    gene_preview = "  ".join(f"`{g}`" for g in condition["genes"][:4])
    if len(condition["genes"]) > 4:
        gene_preview += "  ..."

    with st.expander(f"**{condition['display']}** — {gene_preview}", expanded=True):
        st.markdown(
            '<span class="badge badge-curated">✓ Curated — North Africa</span>',
            unsafe_allow_html=True
        )

        col_genes, col_meta = st.columns([3, 1])
        with col_genes:
            st.markdown('<div class="section-label">Candidate genes</div>', unsafe_allow_html=True)
            chips = " ".join(f'<span class="gene-chip">{g}</span>' for g in condition["genes"])
            st.markdown(chips, unsafe_allow_html=True)
        with col_meta:
            st.markdown(
                f'<span class="badge badge-icd">ICD-10: {condition["icd10"]}</span>',
                unsafe_allow_html=True
            )
            st.write("")
            for i, gene in enumerate(condition["genes"]):
                if i < len(condition["inheritance"]):
                    st.caption(f"`{gene}` — {condition['inheritance'][i]}")

        st.markdown(
            f'<div class="clinical-box">{condition["clinical_note"]}</div>',
            unsafe_allow_html=True
        )

        if condition["key_snps"]:
            st.markdown('<div class="section-label">Known pathogenic variants</div>', unsafe_allow_html=True)
            for rsid, desc in condition["key_snps"].items():
                st.markdown(
                    f'<div class="snp-row"><code>{rsid}</code> &nbsp;—&nbsp; {desc}</div>',
                    unsafe_allow_html=True
                )

        st.markdown('<div class="section-label">Genomic data — pyGeno</div>', unsafe_allow_html=True)
        show_pygeno_section(condition["genes"])


def render_hpo_result(hit: dict):
    with st.expander(f"**{hit['hp_name']}** — `{hit['hp_code']}`", expanded=False):
        st.markdown(
            f'<span class="badge badge-hpo">HPO</span> '
            f'<span style="font-size:12px;color:#3a4868;margin-left:6px;">'
            f'{hit["match_type"]} match</span>',
            unsafe_allow_html=True
        )
        st.markdown('<div class="section-label">Associated genes ({} total)</div>'.format(
            len(hit["genes"])), unsafe_allow_html=True)
        shown = hit["genes"][:20]
        chips = " ".join(f'<span class="gene-chip">{g}</span>' for g in shown)
        if len(hit["genes"]) > 20:
            chips += f'<span style="color:#3a4868;font-size:12px;margin-left:6px;">+{len(hit["genes"])-20} more</span>'
        st.markdown(chips, unsafe_allow_html=True)
        st.markdown('<div class="section-label">Genomic data — pyGeno</div>', unsafe_allow_html=True)
        show_pygeno_section(hit["genes"][:6])


def render_orphanet_result(disease: dict):
    with st.expander(f"**{disease['name']}** — Orphanet {disease['orpha_code']}", expanded=False):
        st.markdown(
            '<span class="badge badge-orphanet">Orphanet</span>',
            unsafe_allow_html=True
        )
        st.markdown('<div class="section-label">Associated genes</div>', unsafe_allow_html=True)
        chips = " ".join(f'<span class="gene-chip">{g}</span>' for g in disease["genes"])
        st.markdown(chips, unsafe_allow_html=True)
        st.markdown(
            f'<a href="https://www.orpha.net/consor/cgi-bin/OC_Exp.php?Expert={disease["orpha_code"]}" '
            f'target="_blank" style="font-size:12px;color:#4a6cf7;text-decoration:none;">'
            f'View on Orphanet →</a>',
            unsafe_allow_html=True
        )
        st.markdown('<div class="section-label">Genomic data — pyGeno</div>', unsafe_allow_html=True)
        show_pygeno_section(disease["genes"][:6])

# ─── SIDEBAR ─────────────────────────────────────────────────────────────────

with st.sidebar:
    st.markdown("## pyGeno Scouter")
    st.caption("Clinical phenotype → genomic coordinates")
    st.divider()

    st.markdown("#### Databases")
    st.success(f"HPO — {len(term_to_id):,} terms")
    st.success(f"HPO genes — {len(hp_to_genes):,} codes")
    st.success(f"Orphanet — {len(orphanet_diseases):,} diseases")
    if pygeno_available():
        st.success(f"pyGeno ({GENOME_BUILD})")
    else:
        st.warning("pyGeno — not configured")
        st.caption("Set PYGENO_PYTHON in app.py")

    st.divider()
    st.markdown("#### Search layers")
    st.caption("**Curated** — North Africa-focused, French/English")
    st.caption("**HPO** — any genetic condition, English")
    st.caption("**Orphanet** — 4,128 rare diseases, French/English")

    st.divider()
    st.markdown("#### Curated conditions")
    for key, entry in PHENOTYPE_DB.items():
        st.caption(f"· {entry['display']}")

    st.divider()
    st.caption("pyGeno · Tariq Daouda, IRIC Montréal")
    st.caption("HPO · JAX &nbsp;|&nbsp; Orphanet · INSERM")
    st.caption(f"Genome: Ensembl {GENOME_BUILD} · Offline")

# ─── HEADER ──────────────────────────────────────────────────────────────────

st.markdown("# 🧬 pyGeno Scouter")
st.markdown(
    '<p style="color:#506080;font-size:14px;margin-top:-8px;margin-bottom:24px;">'
    'Enter a clinical phenotype to identify candidate genes, genomic coordinates, '
    'and known pathogenic variants.</p>',
    unsafe_allow_html=True
)

# ─── SEARCH ──────────────────────────────────────────────────────────────────

col_q, col_btn = st.columns([5, 1])
with col_q:
    query = st.text_input(
        label="Search",
        placeholder="hemolytic anemia · muscle weakness · recurrent fever · G6PD · DMD ...",
        label_visibility="collapsed"
    )
with col_btn:
    st.button("Search", type="primary", use_container_width=True)

# ─── RESULTS ─────────────────────────────────────────────────────────────────

if query and query.strip():
    q = query.strip()

    # Layer 1 — Curated
    curated_matches = search_phenotypes(q)
    if curated_matches:
        st.markdown(
            f'<div class="section-label">Curated results — {len(curated_matches)} match(es)</div>',
            unsafe_allow_html=True
        )
        for condition_key, condition in curated_matches:
            render_curated_result(condition_key, condition)

    # Layer 2 — HPO + Orphanet
    st.markdown('<hr class="divider">', unsafe_allow_html=True)
    with st.expander(
        "Expand search — HPO & Orphanet",
        expanded=(len(curated_matches) == 0)
    ):
        hpo_tab, orphanet_tab = st.tabs(["HPO — Human Phenotype Ontology", "Orphanet"])

        with hpo_tab:
            hpo_hits = search_hpo(q, term_to_id, hp_to_genes, id_to_name, max_results=5)
            if hpo_hits:
                st.caption(f"{len(hpo_hits)} HPO term(s) — HPO is indexed in English")
                for hit in hpo_hits:
                    render_hpo_result(hit)
            else:
                st.info("No HPO terms found. HPO is in English — try: 'muscle weakness', 'hemolytic anemia', 'recurrent fever'")

        with orphanet_tab:
            orphanet_hits = search_orphanet(q, orphanet_diseases, max_results=5)
            if orphanet_hits:
                st.caption(f"{len(orphanet_hits)} Orphanet disease(s) found")
                for disease in orphanet_hits:
                    render_orphanet_result(disease)
            else:
                st.info("No Orphanet diseases found for this query.")

# ─── EMPTY STATE ─────────────────────────────────────────────────────────────

else:
    st.markdown('<hr class="divider">', unsafe_allow_html=True)

    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown('<div class="section-label">Hematology</div>', unsafe_allow_html=True)
        for t in ["hemolytic anemia", "thalassemia", "sickle cell", "G6PD deficiency"]:
            st.caption(f"· {t}")
    with c2:
        st.markdown('<div class="section-label">Neurology</div>', unsafe_allow_html=True)
        for t in ["muscle weakness elevated CK", "peripheral neuropathy", "epilepsy"]:
            st.caption(f"· {t}")
    with c3:
        st.markdown('<div class="section-label">Inflammation / Other</div>', unsafe_allow_html=True)
        for t in ["recurrent fever", "cystic fibrosis", "hypertrophic cardiomyopathy"]:
            st.caption(f"· {t}")

    st.write("")
    st.markdown(
        '<div class="tip-box">'
        'The curated database works in French and English. '
        'HPO search works best in English. '
        'Orphanet works in both. '
        'If a curated result exists it appears first — '
        'use "Expand search" for HPO and Orphanet coverage.'
        '</div>',
        unsafe_allow_html=True
    )

