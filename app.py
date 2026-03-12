"""
app.py — pyGeno Scouter
=======================
Clinical phenotype → genomic coordinates interface.
Curated Morocco-specific database + HPO/Orphanet fallback.

Run:
    streamlit run app.py

Requirements (Python 3.11):
    pip install streamlit pandas

pyGeno queries run via subprocess in Python 3.6 venv.
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

PYGENO_PYTHON = r"C:\pyGeno_env\Scripts\python.exe"
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
    .stApp { background-color: #f7f8fc; }
    .block-container { padding: 2rem 2.5rem; max-width: 1100px; }
    .stTextInput > div > div > input {
        background: #ffffff;
        border: 1.5px solid #d0d5e8;
        border-radius: 10px;
        font-size: 15px;
        padding: 10px 14px;
        color: #1a1f35;
        box-shadow: 0 1px 4px rgba(0,0,0,0.06);
    }
    .stTextInput > div > div > input:focus {
        border-color: #4a6cf7;
        box-shadow: 0 0 0 3px rgba(74,108,247,0.12);
    }
    .gene-chip {
        display: inline-block;
        background: #e8eeff;
        color: #3550c8;
        border: 1px solid #c5d0f5;
        border-radius: 20px;
        padding: 3px 12px;
        margin: 3px;
        font-size: 13px;
        font-weight: 600;
        font-family: monospace;
    }
    .hpo-chip {
        display: inline-block;
        background: #e8f5e8;
        color: #2d6e2d;
        border: 1px solid #b5d9b5;
        border-radius: 6px;
        padding: 2px 10px;
        margin: 2px;
        font-size: 12px;
        font-family: monospace;
    }
    .snp-row {
        background: #f8f9fe;
        border-radius: 6px;
        padding: 6px 10px;
        margin: 4px 0;
        font-size: 13px;
        border-left: 3px solid #4a6cf7;
    }
    .icd-badge {
        background: #f0e8ff;
        color: #6030b0;
        border-radius: 6px;
        padding: 3px 10px;
        font-size: 12px;
        font-family: monospace;
        font-weight: 600;
    }
    .source-badge-curated {
        background: #e8f5e8; color: #2d6e2d;
        border-radius: 6px; padding: 3px 10px;
        font-size: 11px; font-weight: 600;
    }
    .source-badge-hpo {
        background: #e8eeff; color: #3550c8;
        border-radius: 6px; padding: 3px 10px;
        font-size: 11px; font-weight: 600;
    }
    .source-badge-orphanet {
        background: #f5f0ff; color: #6030b0;
        border-radius: 6px; padding: 3px 10px;
        font-size: 11px; font-weight: 600;
    }
    .clinical-box {
        background: #f0f4ff;
        border-left: 4px solid #4a6cf7;
        border-radius: 6px;
        padding: 12px 16px;
        margin: 10px 0;
        color: #2a3060;
        font-size: 14px;
        line-height: 1.6;
    }
    .genomic-pos {
        font-family: monospace;
        font-size: 13px;
        color: #c07020;
        background: #fff8ee;
        border-radius: 6px;
        padding: 4px 10px;
        display: inline-block;
    }
    .section-label {
        font-size: 11px;
        font-weight: 700;
        letter-spacing: 0.08em;
        color: #8090b0;
        text-transform: uppercase;
        margin: 16px 0 6px 0;
    }
    .warn-box {
        background: #fff8ee;
        border-left: 4px solid #f0a050;
        border-radius: 6px;
        padding: 10px 14px;
        color: #805020;
        font-size: 13px;
    }
    h1 { color: #1a1f35 !important; font-weight: 700 !important; }
    h2, h3 { color: #2a3060 !important; }
    .stButton > button { border-radius: 8px; font-weight: 600; }
</style>
""", unsafe_allow_html=True)

# ─── LOAD DATA ───────────────────────────────────────────────────────────────

@st.cache_resource(show_spinner="Chargement des bases de données génomiques...")
def load_all_data():
    term_to_id, id_to_name   = load_hpo_ontology()
    hp_to_genes, gene_to_hps = load_hpo_gene_map()
    orphanet_diseases        = load_orphanet()
    return term_to_id, id_to_name, hp_to_genes, gene_to_hps, orphanet_diseases

term_to_id, id_to_name, hp_to_genes, gene_to_hps, orphanet_diseases = load_all_data()

# ─── PYGENO HELPERS ──────────────────────────────────────────────────────────

@st.cache_data(ttl=3600, show_spinner=False)
def run_pygeno_query(gene_name: str) -> dict:
    if not os.path.exists(PYGENO_PYTHON):
        return {"error": f"pyGeno introuvable: {PYGENO_PYTHON}"}
    try:
        result = subprocess.run(
            [PYGENO_PYTHON, PYGENO_SCRIPT, gene_name, GENOME_BUILD],
            capture_output=True, text=True, timeout=60
        )
        if result.returncode != 0:
            return {"error": result.stderr.strip()}
        return json.loads(result.stdout)
    except subprocess.TimeoutExpired:
        return {"error": "Timeout pyGeno (>60s)"}
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
            '<div class="warn-box">⚠ pyGeno non configuré. '
            'Définir PYGENO_PYTHON dans app.py.</div>',
            unsafe_allow_html=True
        )
        return

    tabs = st.tabs(genes[:6])
    for tab, gene_name in zip(tabs, genes[:6]):
        with tab:
            with st.spinner(f"Interrogation pyGeno — {gene_name}..."):
                data = run_pygeno_query(gene_name)

            if "error" in data:
                st.error(f"pyGeno: {data['error']}")
                continue

            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Chromosome", data["chromosome"])
            c2.metric("Longueur", f"{data['length_bp']:,} bp")
            c3.metric("Transcrits", data["n_transcripts"])
            c4.metric("Isoformes codants", data["n_coding_transcripts"])

            strand = "Forward (+)" if data["strand"] == 1 else "Reverse (−)"
            st.markdown(
                f'<span class="genomic-pos">📍 chr{data["chromosome"]}:'
                f'{data["start"]:,}–{data["end"]:,} · {strand} · {data["ensembl_id"]}</span>',
                unsafe_allow_html=True
            )

            for iso in data.get("isoforms", []):
                with st.expander(
                    f'{iso["transcript_id"]} — {iso["protein_length"]} aa '
                    f'· {iso["coding_exons"]} exons codants'
                ):
                    if iso["exons"]:
                        df = pd.DataFrame(iso["exons"])
                        df.columns = ["Début", "Fin", "CDS (bp)", "5' CDS"]
                        df.index = [f"Exon {i+1}" for i in range(len(df))]
                        df["Début"] = df["Début"].apply(lambda x: f"{x:,}")
                        df["Fin"]   = df["Fin"].apply(lambda x: f"{x:,}")
                        st.dataframe(df, use_container_width=True)

                    if iso.get("protein_seq_60"):
                        st.code(
                            protein_blocks(iso["protein_seq_60"]) +
                            (f"\n... [{iso['protein_length'] - 60} aa supplémentaires]"
                             if iso["protein_length"] > 60 else ""),
                            language=None
                        )

# ─── RESULT RENDERERS ────────────────────────────────────────────────────────

def render_curated_result(condition_key, condition):
    with st.expander(
        f"🇲🇦  {condition['display']}  —  "
        + "  ".join(f"`{g}`" for g in condition["genes"][:4])
        + ("  ..." if len(condition["genes"]) > 4 else ""),
        expanded=True
    ):
        st.markdown(
            '<span class="source-badge-curated">✓ Base Maroc / CHU Ibn Rochd</span>',
            unsafe_allow_html=True
        )

        col_genes, col_meta = st.columns([3, 1])
        with col_genes:
            st.markdown('<div class="section-label">Gènes candidats</div>', unsafe_allow_html=True)
            chips = " ".join(f'<span class="gene-chip">{g}</span>' for g in condition["genes"])
            st.markdown(chips, unsafe_allow_html=True)
        with col_meta:
            st.markdown(
                f'<span class="icd-badge">ICD-10: {condition["icd10"]}</span>',
                unsafe_allow_html=True
            )
            for i, gene in enumerate(condition["genes"]):
                if i < len(condition["inheritance"]):
                    st.caption(f"`{gene}` — {condition['inheritance'][i]}")

        st.markdown(
            f'<div class="clinical-box">📋 <strong>Note clinique:</strong><br>'
            f'{condition["clinical_note"]}</div>',
            unsafe_allow_html=True
        )

        if condition["key_snps"]:
            st.markdown('<div class="section-label">Variants pathogènes connus</div>', unsafe_allow_html=True)
            for rsid, desc in condition["key_snps"].items():
                st.markdown(
                    f'<div class="snp-row"><code>{rsid}</code> — {desc}</div>',
                    unsafe_allow_html=True
                )

        st.markdown('<div class="section-label">Données génomiques (pyGeno)</div>', unsafe_allow_html=True)
        show_pygeno_section(condition["genes"])


def render_hpo_result(hit: dict):
    with st.expander(
        f"🔬  {hit['hp_name']}  —  {hit['hp_code']}",
        expanded=False
    ):
        st.markdown(
            '<span class="source-badge-hpo">HPO — Human Phenotype Ontology</span>',
            unsafe_allow_html=True
        )
        st.markdown(
            f'<span class="hpo-chip">{hit["hp_code"]}</span> '
            f'<span style="font-size:13px;color:#606080;">'
            f'Correspondance: <strong>{hit["match_type"]}</strong></span>',
            unsafe_allow_html=True
        )

        st.markdown('<div class="section-label">Gènes associés ({} total)</div>'.format(
            len(hit["genes"])), unsafe_allow_html=True)
        shown = hit["genes"][:20]
        chips = " ".join(f'<span class="gene-chip">{g}</span>' for g in shown)
        if len(hit["genes"]) > 20:
            chips += f' <span style="color:#8090b0;font-size:12px;">+{len(hit["genes"])-20} autres</span>'
        st.markdown(chips, unsafe_allow_html=True)

        st.markdown('<div class="section-label">Données génomiques (pyGeno) — 6 premiers gènes</div>', unsafe_allow_html=True)
        show_pygeno_section(hit["genes"][:6])


def render_orphanet_result(disease: dict):
    with st.expander(
        f"🏥  {disease['name']}  —  Orphanet: {disease['orpha_code']}",
        expanded=False
    ):
        st.markdown(
            '<span class="source-badge-orphanet">Orphanet</span>',
            unsafe_allow_html=True
        )

        st.markdown('<div class="section-label">Gènes associés</div>', unsafe_allow_html=True)
        chips = " ".join(f'<span class="gene-chip">{g}</span>' for g in disease["genes"])
        st.markdown(chips, unsafe_allow_html=True)

        st.markdown(
            f'<a href="https://www.orpha.net/consor/cgi-bin/OC_Exp.php?Expert={disease["orpha_code"]}" '
            f'target="_blank" style="font-size:13px;color:#4a6cf7;">📖 Voir sur Orphanet →</a>',
            unsafe_allow_html=True
        )

        st.markdown('<div class="section-label">Données génomiques (pyGeno) — 6 premiers gènes</div>', unsafe_allow_html=True)
        show_pygeno_section(disease["genes"][:6])

# ─── SIDEBAR ─────────────────────────────────────────────────────────────────

with st.sidebar:
    st.markdown("## 🧬 pyGeno Scouter")
    st.markdown("**CHU Ibn Rochd · FMPC Casablanca**")
    st.caption("Phénotype clinique → Coordonnées génomiques")
    st.divider()

    st.markdown("#### Bases de données")
    st.success(f"✓ HPO — {len(term_to_id):,} termes")
    st.success(f"✓ HPO genes — {len(hp_to_genes):,} codes")
    st.success(f"✓ Orphanet — {len(orphanet_diseases):,} maladies")
    if pygeno_available():
        st.success(f"✓ pyGeno ({GENOME_BUILD})")
    else:
        st.warning("⚠ pyGeno non configuré")

    st.divider()
    st.markdown("#### Mode de recherche")
    st.caption("**1.** Base Maroc — résultats instantanés, FR")
    st.caption("**2.** HPO — toute pathologie génétique, EN")
    st.caption("**3.** Orphanet — 4128 maladies rares, FR/EN")

    st.divider()
    st.markdown("#### Conditions Maroc")
    for key, entry in PHENOTYPE_DB.items():
        st.caption(f"• {entry['display']}")

    st.divider()
    st.caption("pyGeno — Tariq Daouda, IRIC Montréal")
    st.caption("HPO — JAX · Orphanet — INSERM")
    st.caption("Genome: Ensembl GRCh38 · 100% offline")

# ─── MAIN ─────────────────────────────────────────────────────────────────────

st.title("🧬 pyGeno Scouter")
st.markdown(
    "Saisir un **phénotype clinique** pour identifier les gènes candidats, "
    "les coordonnées génomiques, et les variants pathogènes associés."
)

col_q, col_btn = st.columns([5, 1])
with col_q:
    query = st.text_input(
        label="Recherche",
        placeholder="Ex: anémie hémolytique · muscle weakness · fièvre récurrente · DMD...",
        label_visibility="collapsed"
    )
with col_btn:
    st.button("Rechercher", type="primary", use_container_width=True)

# ─── RESULTS ──────────────────────────────────────────────────────────────────

if query and query.strip():
    q = query.strip()

    # Layer 1 — Curated Morocco DB
    curated_matches = search_phenotypes(q)
    if curated_matches:
        st.markdown(
            f'<div class="section-label">Base Maroc — {len(curated_matches)} résultat(s)</div>',
            unsafe_allow_html=True
        )
        for condition_key, condition in curated_matches:
            render_curated_result(condition_key, condition)

    # Layer 2 — HPO + Orphanet
    st.markdown("---")
    with st.expander(
        "🔬 Élargir la recherche — HPO & Orphanet",
        expanded=(len(curated_matches) == 0)
    ):
        hpo_tab, orphanet_tab = st.tabs(
            ["HPO — Human Phenotype Ontology", "Orphanet"]
        )

        with hpo_tab:
            hpo_hits = search_hpo(q, term_to_id, hp_to_genes, id_to_name, max_results=5)
            if hpo_hits:
                st.caption(f"{len(hpo_hits)} terme(s) HPO correspondant(s) — recherche en anglais recommandée")
                for hit in hpo_hits:
                    render_hpo_result(hit)
            else:
                st.info("Aucun terme HPO trouvé.")
                st.caption("HPO est en anglais. Essayez: 'muscle weakness', 'hemolytic anemia', 'recurrent fever'")

        with orphanet_tab:
            orphanet_hits = search_orphanet(q, orphanet_diseases, max_results=5)
            if orphanet_hits:
                st.caption(f"{len(orphanet_hits)} maladie(s) Orphanet correspondante(s)")
                for disease in orphanet_hits:
                    render_orphanet_result(disease)
            else:
                st.info("Aucune maladie Orphanet trouvée.")

# ─── EMPTY STATE ──────────────────────────────────────────────────────────────

else:
    st.markdown("---")
    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown("**🩸 Hématologie**")
        for t in ["anémie hémolytique", "thalassémie", "drépanocytose", "déficience G6PD"]:
            st.caption(f"• {t}")
    with c2:
        st.markdown("**🧠 Neurologie**")
        for t in ["faiblesse musculaire CK élevée", "neuropathie périphérique", "épilepsie"]:
            st.caption(f"• {t}")
    with c3:
        st.markdown("**🔥 Inflammation / Autre**")
        for t in ["fièvre récurrente", "mucoviscidose", "cardiomyopathie hypertrophique"]:
            st.caption(f"• {t}")

    st.markdown("---")
    st.markdown(
        '<div class="clinical-box">'
        '💡 <strong>Note:</strong> La base Maroc fonctionne en français. '
        'HPO fonctionne mieux en anglais (<em>muscle weakness</em>, '
        '<em>hemolytic anemia</em>). '
        'Orphanet fonctionne dans les deux langues.'
        '</div>',
        unsafe_allow_html=True
    )
