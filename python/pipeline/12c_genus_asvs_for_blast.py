# 12c_genus_asvs_for_blast.py
# The main ASVs of one genus, for BLAST, as Table S2 does for Caproiciproducens
# (12_table_s2_shap_ranking.py): the genus's ASVs that together hold 90% of its
# reads in the modelled Illumina samples (all of its ASVs, before the prevalence
# filter), most reads first. Added 2026-10-14 for Oscillibacter, which the
# Discussion treats as a hypothesis only, so its identity must come from the
# corrected run's primer-free sequences.
#
# No model is fitted, so the output is the same on every machine. The SHAP rank
# of each ASV is read from 12's ranking (data/table_s2/12_asv_shap_ranking.csv);
# an ASV removed by the prevalence filter has none. The Nanopore study has no
# ASVs (its genera come from full-length reads), so it is not covered.
#
# Usage, from the repository root:
#   python python/pipeline/12c_genus_asvs_for_blast.py <Genus> <output folder>
#       writes <Genus>_ASVs_<date>.fasta and <Genus>_asvs.csv (never over
#       existing files);
#   python python/pipeline/12c_genus_asvs_for_blast.py <Genus> <output folder> <BLAST text file>
#       adds the closest relatives from NCBI BLAST's text file (Nucleotide
#       BLAST, database "16S ribosomal RNA sequences (Bacteria and Archaea)",
#       megablast, Download All -> Text), read as 12b reads Table S2's, and
#       writes <Genus>_ASVs_<date>_blast.csv.
# Inputs (run of record in results/primer_trimmed_2026-10-06/): the R/04 matrix
# with every ASV (BIOTWIN_FINAL_ML_MATRIX.csv) and the taxonomy with the
# sequences (ASV_Taxonomy_Master.csv).

import datetime
import glob
import importlib.util
import os
import re
import sys
from pathlib import Path

import pandas as pd

from pipeline_settings import drop_unmodelled

REPO = Path(__file__).resolve().parents[2]
RUN = REPO / "results" / "primer_trimmed_2026-10-06"
ASV_ALL_PATH = RUN / "corrected_2026-10-11" / "BIOTWIN_FINAL_ML_MATRIX.csv"
TAXONOMY_PATH = RUN / "ASV_Taxonomy_Master.csv"
RANKING_PATH = REPO / "data" / "table_s2" / "12_asv_shap_ranking.csv"
GENUS_READ_SHARE = 0.90
# Stage 1's check: no ASV may start with the 341F primer or end with the
# reverse complement of 805R (GACTACHVGGGTATCTAATCC).
PRIMER_341F_START = 'CCTACGGG'
PRIMER_805R_RC_END = re.compile(r'GGATTAGATACCC[CGT][AGT]GTAGTC$')


def load_12b():
    spec = importlib.util.spec_from_file_location("blast_12b", Path(__file__).parent / "12b_table_s2_from_blast.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def genus_asvs(genus_name):
    every = pd.read_csv(ASV_ALL_PATH, low_memory=False)
    for c in ['Paper_ID', 'Sample_ID', 'BIOREACTOR']:
        every[c] = every[c].astype(str).str.strip()
    every = drop_unmodelled(every)
    tax = pd.read_csv(TAXONOMY_PATH).set_index('ASV_ID')
    asvs = [c for c in every.columns if c.startswith('ASV_') and tax['Genus'].get(c) == genus_name]
    if not asvs:
        raise SystemExit(f"No ASV of {ASV_ALL_PATH.name} has the SILVA genus '{genus_name}'.")
    counts = every[asvs].fillna(0)
    reads = counts.sum()
    reads = reads[reads > 0]
    order = sorted(reads.index, key=lambda a: (-reads[a], int(a.split('_')[1])))
    cum = reads[order].cumsum() / reads.sum()
    chosen = order[:int((cum < GENUS_READ_SHARE).sum()) + 1]
    ranking = pd.read_csv(RANKING_PATH).set_index('ASV_ID')

    def rank(a):
        if a not in ranking.index:
            return 'not kept'
        return f"{ranking.loc[a, 'rank']:g}" if ranking.loc[a, 'mean_abs_shap'] > 0 else 'none'

    table = pd.DataFrame({
        'ASV_ID': chosen,
        'pct_of_genus_reads': [100 * reads[a] / reads.sum() for a in chosen],
        'cum_pct_of_genus_reads': [100 * cum[a] for a in chosen],
        'shap_rank': [rank(a) for a in chosen],
        'n_samples': [int((counts[a] > 0).sum()) for a in chosen],
        'n_studies': [every.loc[counts[a] > 0, 'Paper_ID'].nunique() for a in chosen],
        'studies': [', '.join(sorted(every.loc[counts[a] > 0, 'Paper_ID'].unique())) for a in chosen],
        'sequence': tax['Raw_Sequence'].reindex(chosen).values,
    })
    table['length'] = table['sequence'].str.len()
    bad = [a for a, s in zip(table['ASV_ID'], table['sequence'])
           if s.startswith(PRIMER_341F_START) or PRIMER_805R_RC_END.search(s)]
    if bad:
        raise SystemExit(f"Primer sequence found in {bad}; the FASTA must hold primer-free sequences.")
    print(f"{genus_name}: {len(reads)} ASVs with reads in the {len(every)} modelled Illumina samples "
          f"({int(reads.sum())} reads); the {len(chosen)} below hold {table['cum_pct_of_genus_reads'].iloc[-1]:.1f}%.")
    return table


def write_fasta(genus_name, folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    fasta = folder / f"{genus_name}_ASVs_{datetime.date.today():%Y-%m-%d}.fasta"
    csv = folder / f"{genus_name}_asvs.csv"
    for p in (fasta, csv):
        if p.exists():
            raise SystemExit(f"{p} already exists.")
    table = genus_asvs(genus_name)
    print(table.drop(columns='sequence').to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    table.to_csv(csv, index=False)
    with open(fasta, 'w', encoding='utf-8', newline='\n') as f:
        for _, r in table.iterrows():
            f.write(f">{r['ASV_ID']}_{genus_name}_{r['pct_of_genus_reads']:.1f}pct_of_genus_reads\n{r['sequence']}\n")
    print(f"  Saved: {csv}\n  Saved: {fasta}")
    print("Next (by hand): BLAST this FASTA at blast.ncbi.nlm.nih.gov (Nucleotide BLAST, database "
          "'16S ribosomal RNA sequences (Bacteria and Archaea)', program megablast), then Download All -> Text, "
          "and run this script again with that file.")


def add_blast(genus_name, folder, blast_path):
    b = load_12b()
    fastas = sorted(glob.glob(os.path.join(folder, f'{genus_name}_ASVs_????-??-??.fasta')))
    if len(fastas) != 1:
        raise SystemExit(f"Expected one {genus_name}_ASVs_<date>.fasta in {folder}, found {len(fastas)}.")
    fasta = b.read_fasta(fastas[0])
    table = pd.read_csv(os.path.join(folder, f'{genus_name}_asvs.csv'))
    if list(table['ASV_ID']) != list(fasta) or list(table['sequence']) != list(fasta.values()):
        raise SystemExit(f"{fastas[0]} does not match {genus_name}_asvs.csv.")
    out_path = fastas[0].replace('.fasta', '_blast.csv')
    if os.path.exists(out_path):
        raise SystemExit(f"{out_path} already exists.")
    header, queries = b.parse_blast_text(blast_path)
    if header['Program'].split()[0:1] != ['BLASTN'] or '16S_ribosomal_RNA' not in header['Database']:
        raise SystemExit(f"Expected a BLASTN search of the 16S ribosomal RNA database; this file has "
                         f"program '{header['Program']}' and database '{header['Database']}'.")
    by_asv = {}
    for q in queries:
        m = b.ASV_IN_TITLE.match(q['title'])
        if not m or m[1] in by_asv:
            raise SystemExit(f"Query title '{q['title']}' does not name a single new ASV.")
        by_asv[m[1]] = q
    if set(by_asv) != set(fasta) or any(by_asv[a]['length'] != len(fasta[a]) for a in fasta):
        raise SystemExit("The BLAST queries do not match the FASTA's ASVs and lengths.")
    rows = [{**r.drop(['sequence']).to_dict(), **b.summarise(by_asv[r['ASV_ID']]['hits']),
             'n_hits': len(by_asv[r['ASV_ID']]['hits']), 'RID': header['RID']} for _, r in table.iterrows()]
    out = pd.DataFrame(rows)
    out.to_csv(out_path, index=False)
    print(f"BLAST file: {blast_path} (RID {header['RID']}; {header['Database']})")
    print("A hit of 99% or more is the closest known relative, not a species identification.")
    print(f"{'ASV':<9} {'% of genus':>10} {'cumulative':>10} {'SHAP rank':>10} {'studies':>8}  closest relative(s)")
    for _, t in out.iterrows():
        print(f"{t['ASV_ID']:<9} {t['pct_of_genus_reads']:>9.1f}% {t['cum_pct_of_genus_reads']:>9.1f}% "
              f"{t['shap_rank']:>10} {t['n_studies']:>8}  {b.blast_columns(t)}")
    print(f"\n{genus_name} reads in these {len(out)} ASVs, by closest relative:")
    for (rel, hi), g in out.groupby(['closest_relatives', 'at_least_99']):
        print(f"  {rel}, identity {'99% or more' if hi else 'below 99%'}: {len(g)} ASV(s), "
              f"{g['pct_of_genus_reads'].sum():.1f}% of the genus's reads")
    print(f"  Saved: {out_path}")


if __name__ == '__main__':
    if len(sys.argv) == 3:
        write_fasta(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        add_blast(*sys.argv[1:])
    else:
        raise SystemExit("Usage: python 12c_genus_asvs_for_blast.py <Genus> <output folder> [<BLAST text file>]")
