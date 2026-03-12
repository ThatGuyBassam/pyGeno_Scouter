"""
clinvar.py — ClinVar variant fetcher for pyGeno Scouter
========================================================
Queries NCBI E-utilities (free, no API key required) for
pathogenic and likely pathogenic variants by gene name.

Usage:
    from clinvar import fetch_clinvar_variants
    variants = fetch_clinvar_variants("HBB")
"""

import time
import urllib.request
import urllib.parse
import urllib.error
import json
from typing import Optional

# ─── CONSTANTS ───────────────────────────────────────────────────────────────

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Only show these significance labels — everything else is noise
KEEP_SIGNIFICANCE = {
    "pathogenic",
    "likely pathogenic",
    "pathogenic/likely pathogenic",
}

MAX_VARIANTS = 15      # cap per gene — enough to be useful, not overwhelming
REQUEST_DELAY = 0.34   # NCBI asks ≤3 requests/sec without an API key


# ─── HELPERS ─────────────────────────────────────────────────────────────────

def _get(url: str, params: dict) -> Optional[dict]:
    """HTTP GET with NCBI E-utilities, returns parsed JSON or None on error."""
    query = urllib.parse.urlencode(params)
    full_url = f"{url}?{query}"
    try:
        with urllib.request.urlopen(full_url, timeout=10) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception:
        return None


# ─── PUBLIC API ──────────────────────────────────────────────────────────────

def fetch_clinvar_variants(gene_name: str) -> list[dict]:
    """
    Return a list of pathogenic/likely-pathogenic variants for gene_name.

    Each dict has:
        title       str   — variant name / description
        significance str  — clinical significance label
        condition   str   — associated condition(s)
        variation_id str  — ClinVar variation ID (for URL)
        review_status str — review status (stars proxy)
    """
    gene = gene_name.strip().upper()

    # ── Step 1: search ClinVar for this gene + pathogenic filter ─────────────
    search_params = {
        "db": "clinvar",
        "term": (
            f"{gene}[gene] AND ("
            f'"pathogenic"[clinical_significance] OR '
            f'"likely pathogenic"[clinical_significance])'
        ),
        "retmax": MAX_VARIANTS,
        "retmode": "json",
    }
    time.sleep(REQUEST_DELAY)
    search_data = _get(ESEARCH_URL, search_params)

    if not search_data:
        return []

    ids = search_data.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []

    # ── Step 2: fetch summaries for those IDs ────────────────────────────────
    time.sleep(REQUEST_DELAY)
    summary_params = {
        "db": "clinvar",
        "id": ",".join(ids),
        "retmode": "json",
    }
    summary_data = _get(ESUMMARY_URL, summary_params)

    if not summary_data:
        return []

    result_map = summary_data.get("result", {})
    uids = result_map.get("uids", [])

    variants = []
    for uid in uids:
        rec = result_map.get(uid, {})

        # Clinical significance
        sig_obj = rec.get("clinical_significance", {})
        sig = sig_obj.get("description", "").lower().strip()

        if sig not in KEEP_SIGNIFICANCE:
            continue

        # Condition name — flatten list if needed
        trait_set = rec.get("trait_set", [])
        conditions = []
        for trait in trait_set:
            trait_name = trait.get("trait_name", "")
            if trait_name and trait_name.lower() not in ("not provided", "not specified", ""):
                conditions.append(trait_name)
        condition_str = "; ".join(conditions) if conditions else "—"

        # Review status
        review_status = sig_obj.get("review_status", "")

        variants.append({
            "title": rec.get("title", uid),
            "significance": sig.title(),
            "condition": condition_str,
            "variation_id": uid,
            "review_status": review_status,
        })

    return variants


def clinvar_url(variation_id: str) -> str:
    """Direct URL to a ClinVar variation entry."""
    return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/"
