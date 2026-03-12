"""
data_loader.py
==============
Parses and indexes the three local data files:
  - hp.obo              → term name/synonym → HP code
  - phenotype_to_genes  → HP code → gene symbols
  - en_product6.xml     → Orphanet disease → gene symbols

All data stays in memory after first load (cached via Streamlit).
No internet required after initial file download.
"""

import os
import re
import xml.etree.ElementTree as ET

DATA_DIR = r"C:\pyGeno_Scouter\data"

HPO_OBO             = os.path.join(DATA_DIR, "hp.obo")
HPO_PHENOTYPE_GENES = os.path.join(DATA_DIR, "phenotype_to_genes.txt")
ORPHANET_XML        = os.path.join(DATA_DIR, "en_product6.xml")


# ─── HP.OBO PARSER ──────────────────────────────────────────────────────────

def load_hpo_ontology():
    """
    Parse hp.obo and return two dicts:
      term_to_id:  lowercase term name/synonym → HP code
                   e.g. "muscle weakness" → "HP:0001324"
      id_to_name:  HP code → canonical English name
                   e.g. "HP:0001324" → "Muscle weakness"
    """
    term_to_id = {}
    id_to_name = {}

    current_id   = None
    current_name = None
    in_term      = False

    with open(HPO_OBO, encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if line == "[Term]":
                in_term = True
                current_id   = None
                current_name = None
                continue

            if line == "" and in_term:
                in_term = False
                continue

            if not in_term:
                continue

            if line.startswith("id: "):
                current_id = line[4:].strip()

            elif line.startswith("name: ") and current_id:
                current_name = line[6:].strip()
                id_to_name[current_id] = current_name
                term_to_id[current_name.lower()] = current_id

            elif line.startswith("synonym: ") and current_id:
                # synonyms look like: synonym: "Muscle hypotonia" EXACT []
                match = re.search(r'"([^"]+)"', line)
                if match:
                    syn = match.group(1).strip()
                    # Only add if not already mapped to avoid overwriting
                    if syn.lower() not in term_to_id:
                        term_to_id[syn.lower()] = current_id

            elif line.startswith("is_obsolete: true"):
                # Remove obsolete terms
                if current_name and current_name.lower() in term_to_id:
                    del term_to_id[current_name.lower()]

    return term_to_id, id_to_name


# ─── PHENOTYPE_TO_GENES PARSER ───────────────────────────────────────────────

def load_hpo_gene_map():
    """
    Parse phenotype_to_genes.txt and return:
      hp_to_genes: HP code → set of gene symbols
                   e.g. "HP:0001324" → {"DMD", "DYSF", "CAPN3", ...}
      gene_to_hps: gene symbol → set of HP codes
                   e.g. "DMD" → {"HP:0001324", "HP:0003236", ...}
    """
    hp_to_genes  = {}
    gene_to_hps  = {}

    with open(HPO_PHENOTYPE_GENES, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 4:
                continue

            # Columns: hpo_id, hpo_name, ncbi_gene_id, gene_symbol, disease_id
            hp_code     = parts[0].strip()
            gene_symbol = parts[3].strip()

            if not hp_code.startswith("HP:") or not gene_symbol:
                continue

            if hp_code not in hp_to_genes:
                hp_to_genes[hp_code] = set()
            hp_to_genes[hp_code].add(gene_symbol)

            if gene_symbol not in gene_to_hps:
                gene_to_hps[gene_symbol] = set()
            gene_to_hps[gene_symbol].add(hp_code)

    return hp_to_genes, gene_to_hps


# ─── ORPHANET XML PARSER ─────────────────────────────────────────────────────

def load_orphanet():
    """
    Parse en_product6.xml and return:
      orphanet_diseases: list of dicts, each with:
        - orpha_code:  str  e.g. "166024"
        - name:        str  e.g. "Duchenne muscular dystrophy"
        - genes:       list of gene symbols
        - name_lower:  str  lowercased name for search
    """
    tree = ET.parse(ORPHANET_XML)
    root = tree.getroot()

    diseases = []

    for disorder in root.iter("Disorder"):
        orpha_code_el = disorder.find("OrphaCode")
        name_el       = disorder.find("Name")

        if orpha_code_el is None or name_el is None:
            continue

        orpha_code = orpha_code_el.text.strip()
        name       = name_el.text.strip()

        # Collect all associated gene symbols
        genes = []
        for gene_el in disorder.iter("Gene"):
            symbol_el = gene_el.find("Symbol")
            if symbol_el is not None and symbol_el.text:
                sym = symbol_el.text.strip()
                if sym and sym not in genes:
                    genes.append(sym)

        if genes:  # Only include diseases with known genes
            diseases.append({
                "orpha_code":  orpha_code,
                "name":        name,
                "name_lower":  name.lower(),
                "genes":       genes,
            })

    return diseases


# ─── COMBINED SEARCH ────────────────────────────────────────────────────────

def search_hpo(query: str, term_to_id: dict, hp_to_genes: dict, id_to_name: dict,
               max_results: int = 8) -> list:
    """
    Search HPO ontology by free text query.
    Returns list of dicts ranked by match quality:
      - hp_code:    HP:XXXXXXX
      - hp_name:    canonical term name
      - genes:      list of associated gene symbols
      - match_type: "exact", "startswith", or "partial"
    """
    query_lower = query.lower().strip()
    results     = []
    seen_hp     = set()

    # Pass 1 — exact match
    if query_lower in term_to_id:
        hp_code = term_to_id[query_lower]
        if hp_code not in seen_hp:
            genes = sorted(hp_to_genes.get(hp_code, set()))
            if genes:
                results.append({
                    "hp_code":    hp_code,
                    "hp_name":    id_to_name.get(hp_code, query),
                    "genes":      genes,
                    "match_type": "exact",
                })
                seen_hp.add(hp_code)

    # Pass 2 — startswith
    for term, hp_code in term_to_id.items():
        if hp_code in seen_hp:
            continue
        if term.startswith(query_lower):
            genes = sorted(hp_to_genes.get(hp_code, set()))
            if genes:
                results.append({
                    "hp_code":    hp_code,
                    "hp_name":    id_to_name.get(hp_code, term),
                    "genes":      genes,
                    "match_type": "startswith",
                })
                seen_hp.add(hp_code)
        if len(results) >= max_results:
            break

    # Pass 3 — partial match (any word in query appears in term)
    if len(results) < max_results:
        words = [w for w in query_lower.split() if len(w) >= 4]
        for term, hp_code in term_to_id.items():
            if hp_code in seen_hp:
                continue
            if any(w in term for w in words):
                genes = sorted(hp_to_genes.get(hp_code, set()))
                if genes:
                    results.append({
                        "hp_code":    hp_code,
                        "hp_name":    id_to_name.get(hp_code, term),
                        "genes":      genes,
                        "match_type": "partial",
                    })
                    seen_hp.add(hp_code)
            if len(results) >= max_results:
                break

    return results[:max_results]


def search_orphanet(query: str, orphanet_diseases: list,
                    max_results: int = 8) -> list:
    """
    Search Orphanet diseases by name.
    Returns list of matching disease dicts ranked by match quality.
    """
    query_lower = query.lower().strip()
    exact       = []
    starts      = []
    partial     = []

    for disease in orphanet_diseases:
        name_lower = disease["name_lower"]
        if name_lower == query_lower:
            exact.append(disease)
        elif name_lower.startswith(query_lower):
            starts.append(disease)
        elif query_lower in name_lower:
            partial.append(disease)
        else:
            # Word-level match
            words = [w for w in query_lower.split() if len(w) >= 4]
            if words and any(w in name_lower for w in words):
                partial.append(disease)

    ranked = exact + starts + partial
    return ranked[:max_results]
