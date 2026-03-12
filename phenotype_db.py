"""
phenotype_db.py
Phenotype-to-Gene mapping database for the pyGeno Scouter.

Maps clinical presentations (symptoms, lab findings, clinical context) to
candidate genes and known pathogenic variants. Focused on conditions relevant
to CHU Ibn Rochd / Morocco patient population.

Each entry contains:
  - genes: list of candidate gene symbols
  - omim: OMIM disease IDs for reference
  - inheritance: AD=autosomal dominant, AR=autosomal recessive, XL=X-linked
  - key_snps: dict of rsID → clinical description
  - clinical_note: what a clinician should know
  - icd10: ICD-10 code for the condition
"""

PHENOTYPE_DB = {

    # ── HEMATOLOGY ──────────────────────────────────────────────────────────

    "anémie hémolytique": {
        "display": "Anémie Hémolytique",
        "genes": ["G6PD", "HBB", "PKLR", "SPTA1"],
        "omim": ["300908", "603903", "266200", "182870"],
        "inheritance": ["XL", "AR", "AR", "AD/AR"],
        "key_snps": {
            "rs1050828": "G6PD A- (Val68Met) — déficience G6PD, fréquent au Maroc",
            "rs5030868": "G6PD Méditerranéenne (Ser218Phe) — forme sévère, Maghreb",
            "rs334":     "HbS (Glu6Val) — drépanocytose",
        },
        "clinical_note": (
            "Chez un patient marocain avec anémie hémolytique, éliminer en priorité "
            "la déficience G6PD (fèves, médicaments déclencheurs) et la bêta-thalassémie. "
            "Demander: hémogramme, réticulocytes, LDH, bilirubine indirecte, test de Coombs."
        ),
        "icd10": "D55–D58",
    },

    "drépanocytose": {
        "display": "Drépanocytose / Syndrome Falciforme",
        "genes": ["HBB"],
        "omim": ["603903"],
        "inheritance": ["AR"],
        "key_snps": {
            "rs334":      "HbS — Glu6Val — mutation causale principale",
            "rs11549407": "HbC — Glu6Lys — double hétérozygotie HbSC possible",
        },
        "clinical_note": (
            "Mutation ponctuelle c.20A>T dans HBB. Homozygotes HbSS = drépanocytose classique. "
            "Hétérozygotes HbAS = trait drépanocytaire (asymptomatique). "
            "Crises vaso-occlusives, AVC, splénomégalie fonctionnelle."
        ),
        "icd10": "D57",
    },

    "thalassémie": {
        "display": "Bêta-Thalassémie",
        "genes": ["HBB", "HBD"],
        "omim": ["613985"],
        "inheritance": ["AR", "AR"],
        "key_snps": {
            "rs35004220": "IVS-I-110 G>A — mutation la plus fréquente au Maroc/Méditerranée",
            "rs35979231": "CD39 C>T — mutation stop, fréquente en Afrique du Nord",
        },
        "clinical_note": (
            "Bêta-thalassémie majeure (Cooley): anémie profonde, hépato-splénomégalie, "
            "déformations osseuses, dépendance transfusionnelle. "
            "Dépistage: VGM bas, hypochromie, électrophorèse de l'hémoglobine (HbA2 élevée)."
        ),
        "icd10": "D56.1",
    },

    "déficience G6PD": {
        "display": "Déficience en G6PD",
        "genes": ["G6PD"],
        "omim": ["300908"],
        "inheritance": ["XL"],
        "key_snps": {
            "rs1050828": "G6PD A- (Val68Met) — activité résiduelle ~12%",
            "rs5030868": "G6PD Méditerranéenne — activité quasi nulle, forme la plus sévère",
            "rs2230037": "G6PD A (Asn126Asp) — activité normale, polymorphisme",
        },
        "clinical_note": (
            "Touche ~5-10% des hommes marocains. Asymptomatique au repos. "
            "Crises déclenchées par: fèves (favisme), primaquine, dapsone, infections. "
            "Éviter les oxydants. Ne pas transfuser en phase aiguë si stable."
        ),
        "icd10": "D55.0",
    },

    # ── NEUROLOGIE / NEUROMUSCULAIRE ────────────────────────────────────────

    "faiblesse musculaire CK élevée": {
        "display": "Faiblesse Musculaire + CK Élevée (Myopathie)",
        "genes": ["DMD", "DYSF", "CAPN3", "ANO5", "SGCA"],
        "omim": ["310200", "253601", "253600", "611307", "600119"],
        "inheritance": ["XL", "AR", "AR", "AR", "AR"],
        "key_snps": {
            "rs80338799": "DMD — délétion exons 45-55 (hot spot Duchenne)",
        },
        "clinical_note": (
            "CK > 10x normale oriente vers dystrophie musculaire. "
            "Duchenne (DMD, garçons <5 ans): faiblesse proximale, signe de Gowers. "
            "Dysferlinopathie: début adulte, atteinte distale ou ceinture. "
            "Bilan: IRM musculaire, biopsie, panel génétique NGS."
        ),
        "icd10": "G71.0",
    },

    "neuropathie périphérique": {
        "display": "Neuropathie Périphérique Héréditaire",
        "genes": ["PMP22", "MPZ", "GJB1", "MFN2"],
        "omim": ["118220", "118200", "302800", "607736"],
        "inheritance": ["AD", "AD", "XL", "AD"],
        "key_snps": {
            "rs104893896": "PMP22 — duplication 17p12 — CMT1A (Charcot-Marie-Tooth type 1A)",
        },
        "clinical_note": (
            "CMT1A: neuropathie démyélinisante, pieds creux, steppage. "
            "EMG: ralentissement des vitesses de conduction. "
            "Duplication PMP22 sur chr17 détectable par MLPA ou CGH-array. "
            "Pas de traitement curatif, kiné et orthèses."
        ),
        "icd10": "G60.0",
    },

    "épilepsie": {
        "display": "Épilepsie Génétique",
        "genes": ["SCN1A", "KCNQ2", "SCN2A", "CDKL5", "PCDH19"],
        "omim": ["607208", "613720", "613721", "300203", "300088"],
        "inheritance": ["AD", "AD", "AD", "XL", "XL"],
        "key_snps": {
            "rs121917959": "SCN1A — syndrome de Dravet (épilepsie sévère du nourrisson)",
        },
        "clinical_note": (
            "Syndrome de Dravet: convulsions fébriles prolongées avant 1 an, "
            "puis épilepsie pharmaco-résistante. Éviter lamotrigine et carbamazépine (aggravent). "
            "Traitement: valproate + clobazam ± stiripentol."
        ),
        "icd10": "G40",
    },

    # ── IMMUNOLOGIE / INFLAMMATION ──────────────────────────────────────────

    "fièvre récurrente": {
        "display": "Fièvre Récurrente / Fièvre Méditerranéenne Familiale",
        "genes": ["MEFV", "MVK", "TNFRSF1A", "NLRP3"],
        "omim": ["249100", "260920", "142680", "191900"],
        "inheritance": ["AR", "AR", "AD", "AD"],
        "key_snps": {
            "rs61752717": "MEFV M694V — mutation FMF la plus sévère, fréquente au Maghreb",
            "rs28940580": "MEFV M680I — FMF forme modérée",
            "rs61752718": "MEFV V726A — FMF forme légère",
            "rs11466023": "MEFV E148Q — variant de signification incertaine",
        },
        "clinical_note": (
            "FMF: crises de 12-72h de fièvre + sérosite (péritonite, pleurite, arthrite). "
            "Fréquence élevée dans les populations du Maghreb, Moyen-Orient, Turquie. "
            "Risque d'amylose AA si non traitée. Traitement: colchicine à vie."
        ),
        "icd10": "E85.0 / M04.1",
    },

    "déficit immunitaire": {
        "display": "Déficit Immunitaire Primitif",
        "genes": ["BTK", "RAG1", "RAG2", "ADA", "IL2RG"],
        "omim": ["300300", "601457", "601457", "102700", "300400"],
        "inheritance": ["XL", "AR", "AR", "AR", "XL"],
        "key_snps": {
            "rs80338705": "BTK — agammaglobulinémie de Bruton (garçons, infections bactériennes récurrentes)",
        },
        "clinical_note": (
            "Agammaglobulinémie de Bruton: garçons avec infections bactériennes récurrentes "
            "après 6 mois (disparition des anticorps maternels). "
            "Bilan: dosage immunoglobulines, sous-populations lymphocytaires, BTK en cytométrie. "
            "Traitement: immunoglobulines IV ou SC à vie."
        ),
        "icd10": "D83.9",
    },

    # ── MÉTABOLISME ─────────────────────────────────────────────────────────

    "phénylcétonurie": {
        "display": "Phénylcétonurie (PKU)",
        "genes": ["PAH"],
        "omim": ["261600"],
        "inheritance": ["AR"],
        "key_snps": {
            "rs5030850": "PAH IVS12+1G>A — mutation fréquente en Afrique du Nord",
        },
        "clinical_note": (
            "Déficit en phénylalanine hydroxylase → accumulation de phénylalanine → "
            "retard mental si non traité. Dépistage néonatal: Guthrie test. "
            "Régime hypoprotidique strict + phénylalanine contrôlée."
        ),
        "icd10": "E70.0",
    },

    "diabète monogénique": {
        "display": "Diabète Monogénique (MODY)",
        "genes": ["GCK", "HNF1A", "HNF4A", "HNF1B", "INS"],
        "omim": ["125851", "600496", "125850", "189907", "222100"],
        "inheritance": ["AD", "AD", "AD", "AD", "AR"],
        "key_snps": {},
        "clinical_note": (
            "MODY suspecté: diabète chez jeune <35 ans, non obèse, ATCD familiaux, "
            "pas d'anticorps anti-ilots. GCK-MODY: hyperglycémie légère stable, souvent "
            "découverte fortuite. HNF1A-MODY: sensible aux sulfamides. "
            "Important de distinguer du T1DM et T2DM car traitement différent."
        ),
        "icd10": "E13",
    },

    # ── FIBROSE KYSTIQUE ────────────────────────────────────────────────────

    "mucoviscidose": {
        "display": "Mucoviscidose / Fibrose Kystique",
        "genes": ["CFTR"],
        "omim": ["219700"],
        "inheritance": ["AR"],
        "key_snps": {
            "rs113993960": "F508del — mutation la plus fréquente (mondiale)",
            "rs75961395": "G542X — fréquente en Méditerranée/Afrique du Nord",
            "rs80034486": "N1303K — prévalence accrue chez patients marocains",
        },
        "clinical_note": (
            "Spectre particulier au Maghreb: F508del moins dominant que dans les populations "
            "européennes, mutations locales importantes. Dépistage: test de la sueur (Cl- >60 mEq/L). "
            "Manifestations: bronchopneumopathies, insuffisance pancréatique exocrine, infertilité masculine."
        ),
        "icd10": "E84",
    },

    # ── CARDIOLOGIE ─────────────────────────────────────────────────────────

    "cardiomyopathie hypertrophique": {
        "display": "Cardiomyopathie Hypertrophique",
        "genes": ["MYH7", "MYBPC3", "TNNT2", "TNNI3"],
        "omim": ["192600", "600958", "115195", "611880"],
        "inheritance": ["AD", "AD", "AD", "AD"],
        "key_snps": {
            "rs121913632": "MYH7 Arg403Gln — phénotype sévère, risque MSC élevé",
            "rs180204588": "MYBPC3 — mutation la plus fréquente toutes ethnies confondues",
        },
        "clinical_note": (
            "CMH: hypertrophie VG asymétrique, non dilatée, non ischémique. "
            "Cause principale de mort subite cardiaque chez le jeune sportif. "
            "ECG: HVG, ondes Q profondes en latéral. Écho: épaisseur septum >15mm. "
            "Éviter sports de compétition. Discussion défibrillateur implantable."
        ),
        "icd10": "I42.1",
    },

    "syndrome du QT long": {
        "display": "Syndrome du QT Long",
        "genes": ["KCNQ1", "KCNH2", "SCN5A"],
        "omim": ["192500", "613688", "603830"],
        "inheritance": ["AD", "AD", "AD"],
        "key_snps": {
            "rs199472709": "KCNQ1 — LQT1, déclenché par effort/nage",
            "rs199473130": "KCNH2 — LQT2, déclenché par bruits soudains",
        },
        "clinical_note": (
            "QTc >480ms chez femme, >470ms chez homme = prolongé. "
            "LQT1: éviter nage, bétabloquants efficaces. "
            "LQT2: éviter bruits soudains, hypokaliémie. "
            "LQT3: mexilétine. Risque torsades de pointes → FV → mort subite."
        ),
        "icd10": "I45.81",
    },
}

# Synonym mapping — maps French clinical terms to database keys
SYNONYMS = {
    "anémie": "anémie hémolytique",
    "hemolytique": "anémie hémolytique",
    "hémolytique": "anémie hémolytique",
    "sicklecell": "drépanocytose",
    "sickle": "drépanocytose",
    "drepanocytose": "drépanocytose",
    "falciforme": "drépanocytose",
    "hbs": "drépanocytose",
    "thalassemie": "thalassémie",
    "thalassemia": "thalassémie",
    "g6pd": "déficience G6PD",
    "fave": "déficience G6PD",
    "favisme": "déficience G6PD",
    "myopathie": "faiblesse musculaire CK élevée",
    "ck": "faiblesse musculaire CK élevée",
    "creatine kinase": "faiblesse musculaire CK élevée",
    "duchenne": "faiblesse musculaire CK élevée",
    "dmd": "faiblesse musculaire CK élevée",
    "cmt": "neuropathie périphérique",
    "charcot": "neuropathie périphérique",
    "pied creux": "neuropathie périphérique",
    "convulsion": "épilepsie",
    "seizure": "épilepsie",
    "dravet": "épilepsie",
    "fmf": "fièvre récurrente",
    "mefv": "fièvre récurrente",
    "peritonite": "fièvre récurrente",
    "peritonitis": "fièvre récurrente",
    "bruton": "déficit immunitaire",
    "agammaglobulinemie": "déficit immunitaire",
    "pku": "phénylcétonurie",
    "phenylcetonurie": "phénylcétonurie",
    "mody": "diabète monogénique",
    "diabete monogenique": "diabète monogénique",
    "mucoviscidose": "mucoviscidose",
    "cftr": "mucoviscidose",
    "cystic fibrosis": "mucoviscidose",
    "cmh": "cardiomyopathie hypertrophique",
    "cardiomyopathie": "cardiomyopathie hypertrophique",
    "hcm": "cardiomyopathie hypertrophique",
    "qt long": "syndrome du QT long",
    "lqt": "syndrome du QT long",
    "torsade": "syndrome du QT long",
}


def search_phenotypes(query: str) -> list:
    """
    Search phenotype DB by symptom keywords.
    Returns list of matching (key, entry) tuples, ranked by relevance.
    """
    query_lower = query.lower().strip()
    matches = {}

    # Direct key match
    for key in PHENOTYPE_DB:
        if query_lower in key.lower():
            matches[key] = 3  # highest weight

    # Synonym match
    for synonym, target_key in SYNONYMS.items():
        if synonym in query_lower or query_lower in synonym:
            matches[target_key] = max(matches.get(target_key, 0), 2)

    # Word-level partial match against display name and clinical note
    words = query_lower.split()
    for key, entry in PHENOTYPE_DB.items():
        searchable = (
            entry["display"].lower() + " " +
            entry["clinical_note"].lower() + " " +
            " ".join(entry["genes"]).lower()
        )
        for word in words:
            if len(word) >= 3 and word in searchable:
                matches[key] = max(matches.get(key, 0), 1)

    # Sort by weight descending
    ranked = sorted(matches.items(), key=lambda x: x[1], reverse=True)
    return [(k, PHENOTYPE_DB[k]) for k, _ in ranked]
