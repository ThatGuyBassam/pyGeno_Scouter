"""
pygeno_query.py — pyGeno query backend
======================================
Runs inside the Python 3.6 pyGeno virtual environment.
Called by the Streamlit app as a subprocess. Accepts a gene name,
queries pyGeno, and returns JSON to stdout.

Usage:
    C:\pyGeno_env\Scripts\python.exe pygeno_query.py HBB GRCh38.78
    C:\pyGeno_env\Scripts\python.exe pygeno_query.py G6PD GRCh38.78

Output: JSON to stdout, errors to stderr.
"""

import sys
import json

def query_gene(gene_name, genome_build):
    try:
        from pyGeno.Genome import Genome
        from pyGeno.Gene import Gene
        from pyGeno.Transcript import Transcript
        from pyGeno.Protein import Protein
        from pyGeno.Exon import Exon
    except ImportError as e:
        return {"error": f"pyGeno not installed in this environment: {e}"}

    try:
        g = Genome(name=genome_build)
    except Exception as e:
        return {"error": f"Genome '{genome_build}' not found. Run setup_genome.py first. ({e})"}

    results = g.get(Gene, name=gene_name)
    if not results:
        return {"error": f"Gene '{gene_name}' not found in {genome_build}"}

    gene = results[0]

    all_transcripts = gene.get(Transcript)
    coding_transcripts = [t for t in all_transcripts if t.biotype == "protein_coding"]

    isoforms = []
    for transcript in coding_transcripts[:3]:  # limit to 3 isoforms
        exons = transcript.get(Exon)
        coding_exons = [e for e in exons if e.CDS]

        proteins = transcript.get(Protein)
        protein_id = proteins[0].id if proteins else None
        protein_seq = proteins[0].sequence if proteins else None
        protein_len = len(protein_seq) if protein_seq else 0

        total_cds = sum(len(e.CDS) for e in coding_exons if e.CDS)

        exon_list = []
        for exon in coding_exons:
            exon_list.append({
                "start":      exon.start,
                "end":        exon.end,
                "cds_length": len(exon.CDS) if exon.CDS else 0,
                "cds_5prime": exon.CDS[:15] if exon.CDS else "",
            })

        isoforms.append({
            "transcript_id":   transcript.id,
            "biotype":         transcript.biotype,
            "total_exons":     len(exons),
            "coding_exons":    len(coding_exons),
            "total_cds_bp":    total_cds,
            "expected_aa":     total_cds // 3,
            "protein_id":      protein_id,
            "protein_length":  protein_len,
            "protein_seq_60":  protein_seq[:60] if protein_seq else None,
            "exons":           exon_list,
        })

    return {
        "gene_name":    gene.name,
        "ensembl_id":   gene.id,
        "biotype":      gene.biotype,
        "chromosome":   str(gene.chromosome.number),
        "start":        gene.start,
        "end":          gene.end,
        "length_bp":    gene.end - gene.start,
        "strand":       gene.strand,
        "genome_build": genome_build,
        "n_transcripts":        len(all_transcripts),
        "n_coding_transcripts": len(coding_transcripts),
        "isoforms":             isoforms,
    }


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(json.dumps({"error": "Usage: pygeno_query.py GENE_NAME GENOME_BUILD"}))
        sys.exit(1)

    gene_name    = sys.argv[1]
    genome_build = sys.argv[2]

    result = query_gene(gene_name, genome_build)
    print(json.dumps(result, indent=2))
