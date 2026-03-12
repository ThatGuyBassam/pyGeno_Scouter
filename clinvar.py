"""
clinvar.py — ClinVar variant fetcher for pyGeno Scouter
========================================================
Queries NCBI E-utilities (free, no API key required) for
pathogenic and likely pathogenic variants by gene name.

Usage:
    from clinvar import fetch_clinvar_variants, clinvar_url
    variants = fetch_clinvar_variants("HBB")
"""

import time
import urllib.request
import urllib.parse
import json
from typing import Optional

# ─── CONSTANTS ───────────────────────────────────────────────────────────────

ESEARCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

# Significance strings to keep — checked against germline_classification.description
KEEP_SIGNIFICANCE = {
    "pathogenic",
    "likely pathogenic",
    "pathogenic/likely pathogenic",
    "pathogenic, other",
    "pathogenic; other",
    "likely pathogenic, other",
}

MAX_FETCH     = 40     # fetch more, filter locally
REQUEST_DELAY = 0.4    # NCBI: <=3 req/sec without API key


# ─── HELPERS ─────────────────────────────────────────────────────────────────

def _get(url: str, params: dict) -> Optional[dict]:
    params["retmode"] = "json"
    full_url = url + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(
        full_url,
        headers={"User-Agent": "pyGenoScouter/1.0 (research tool)"}
    )
    try:
        with urllib.request.urlopen(req, timeout=12) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception:
        return None


# ─── PUBLIC API ──────────────────────────────────────────────────────────────

def fetch_clinvar_variants(gene_name: str) -> list:
    """
    Return pathogenic / likely-pathogenic ClinVar variants for gene_name.

    Each dict:
        title        str  — variant name e.g. NM_000518.5(HBB):c.20A>T
        significance str  — e.g. "Pathogenic"
        condition    str  — associated condition(s)
        variation_id str  — ClinVar variation ID
        review_status str — review status string
    """
    gene = gene_name.strip().upper()

    # ── Step 1: search ───────────────────────────────────────────────────────
    time.sleep(REQUEST_DELAY)
    search_data = _get(ESEARCH_URL, {
        "db":     "clinvar",
        "term":   f"{gene}[gene] AND (pathogenic[clinsig] OR likely pathogenic[clinsig])",
        "retmax": MAX_FETCH,
    })

    if not search_data:
        return []

    ids = search_data.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []

    # ── Step 2: fetch summaries ───────────────────────────────────────────────
    time.sleep(REQUEST_DELAY)
    summary_data = _get(ESUMMARY_URL, {
        "db": "clinvar",
        "id": ",".join(ids),
    })

    if not summary_data:
        return []

    result_map = summary_data.get("result", {})
    uids       = result_map.get("uids", [])

    variants = []
    for uid in uids:
        rec = result_map.get(uid, {})

        # ── Significance — now under germline_classification ──────────────
        germline = rec.get("germline_classification", {})
        sig_raw  = germline.get("description", "").strip()
        sig      = sig_raw.lower()

        if not any(keep in sig for keep in KEEP_SIGNIFICANCE):
            continue

        # ── Condition — from germline_classification.trait_set ────────────
        trait_set  = germline.get("trait_set", [])
        conditions = []
        for trait in trait_set:
            name = trait.get("trait_name", "").strip()
            if name and name.lower() not in ("not provided", "not specified", ""):
                conditions.append(name)
        condition_str = "; ".join(conditions) if conditions else "—"

        variants.append({
            "title":         rec.get("title", uid),
            "significance":  sig_raw.title(),
            "condition":     condition_str,
            "variation_id":  uid,
            "review_status": germline.get("review_status", ""),
        })

    return variants


def clinvar_url(variation_id: str) -> str:
    return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/"

