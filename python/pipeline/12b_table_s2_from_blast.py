# 12b_table_s2_from_blast.py
# Stage 7 of the primer-trimmed rerun (handoff: "Rebuild Table S2"), step 4:
# build Table S2 directly from NCBI BLAST's text file, so no name or number is
# ever copied by hand.
#
# Usage, from the repository root:
#   python python/pipeline/12b_table_s2_from_blast.py <folder> <BLAST text file>
# <folder> holds the outputs of 12_table_s2_shap_ranking.py
# (12_table_s2_asvs.csv and Table_S2_ASVs_<date>.fasta). The BLAST file is the
# web BLAST result for that FASTA (Nucleotide BLAST, database "16S ribosomal RNA
# sequences (Bacteria and Archaea)", megablast), saved with Download All -> Text
# (<RID>-Alignment.txt).
# Output: TableS2_<date>.csv in <folder>, never over an existing file, with one
# row per ASV in the FASTA's order: the top 20 by mean |SHAP|, then the other
# Caproiciproducens ASVs that together hold 90% of the genus's reads
# (included_as says which list an ASV is in).
#
# For each ASV:
#   - the closest relative(s): the top hit, and every other species whose hit
#     has the same max score, identity and query cover (ties are all listed);
#   - their query cover, identity and accession;
#   - whether the identity is at least 99%. Such a hit is the closest known
#     relative, not a species identification;
#   - the next best species outside the tie, with its identity;
#   - from 12: its SHAP rank, the number of samples and studies it is found
#     in, and its share of its genus's reads.
# For the Caproiciproducens ASVs, it also adds up the genus's reads by closest
# relative, so the strain-level evidence covers the whole genus.
#
# Checks: the file is a BLASTN search of the 16S ribosomal RNA database, its
# queries are exactly the FASTA's ASVs, with the same lengths, and the FASTA is
# the one 12 wrote with 12_table_s2_asvs.csv.

import glob
import os
import re
import sys

import pandas as pd

HIT_ROW = re.compile(r'(?P<taxid>\S+)\s+(?P<max>[\d.]+)\s+(?P<total>[\d.]+)\s+(?P<cover>\d+)%\s+'
                     r'(?P<evalue>\S+)\s+(?P<ident>[\d.]+)\s+(?P<acc_len>\d+)\s+(?P<acc>\S+)\s*$')
ASV_IN_TITLE = re.compile(r'(ASV_\d+)(?!\d)')  # titles start with the ASV, e.g. ASV_143_rank_4_Caproiciproducens


def species_of(description):
    """Genus and species from a hit's full description, e.g. 'Caproicibacterium
    lactatifermentans strain LBM19010 16S ribosomal RNA, partial sequence'."""
    words = description.split()
    n = 3 if words[0] == 'Candidatus' else 2
    if len(words) > n + 1 and words[n] in ('subsp.', 'pv.'):
        n += 2
    return ' '.join(words[:n])


def parse_blast_text(path):
    """Header fields, and for each query its title, length and hits in BLAST's order."""
    with open(path, encoding='utf-8', errors='replace') as f:
        text = f.read()
    header = {k: (re.search(rf'^{k}: (.*?)\s*$', text, re.M) or [None, ''])[1]
              for k in ('RID', 'Program', 'Database')}
    blocks = re.split(r'^(?=Query #\d+: )', text, flags=re.M)[1:]
    if not blocks:
        raise SystemExit(f"{path} has no 'Query #' sections. Use BLAST's Download All -> Text, "
                         "not the hit table.")
    queries = []
    for block in blocks:
        q = re.match(r'Query #(\d+): (.*?) Query ID: \S+ Length: (\d+)', block)
        if not q:
            raise SystemExit(f"Cannot read this query line: {block.splitlines()[0]}")
        # Full descriptions, from the alignments (the hit table cuts them short)
        full = {acc: desc for desc, acc in re.findall(r'^>(.*?)\s*\nSequence ID: (\S+)', block, re.M)}
        hits = []
        if 'Sequences producing significant alignments:' in block:
            table = block.split('Sequences producing significant alignments:', 1)[1].split('Alignments:', 1)[0]
            for line in table.splitlines()[3:]:
                if not line.strip():
                    break
                m = HIT_ROW.search(line)
                if not m:
                    raise SystemExit(f"Cannot read this hit row: {line}")
                acc = m['acc']
                if acc not in full:
                    raise SystemExit(f"Hit {acc} has no alignment section in {path}.")
                hits.append({'accession': acc, 'description': full[acc], 'species': species_of(full[acc]),
                             'taxid': m['taxid'], 'max_score': float(m['max']), 'cover': int(m['cover']),
                             'evalue': m['evalue'], 'identity': float(m['ident'])})
        queries.append({'n': int(q[1]), 'title': q[2].strip(), 'length': int(q[3]), 'hits': hits})
    return header, queries


def summarise(hits):
    """The tied top species (all listed) and the next best species outside the tie."""
    if not hits:
        return {'closest_relatives': 'no significant hit', 'query_cover_pct': None, 'identity_pct': None,
                'accessions': '', 'n_species_tied': 0, 'at_least_99': False,
                'next_species': '', 'next_identity_pct': None}
    top = hits[0]
    tied, seen = [], set()
    for h in hits:
        if (h['max_score'], h['identity'], h['cover']) == (top['max_score'], top['identity'], top['cover']) \
                and h['species'] not in seen:
            tied.append(h)
            seen.add(h['species'])
    nxt = next((h for h in hits if h['species'] not in seen), None)
    return {'closest_relatives': ' / '.join(h['species'] for h in tied),
            'query_cover_pct': top['cover'], 'identity_pct': top['identity'],
            'accessions': ' / '.join(h['accession'] for h in tied),
            'n_species_tied': len(tied), 'at_least_99': top['identity'] >= 99,
            'next_species': nxt['species'] if nxt else '',
            'next_identity_pct': nxt['identity'] if nxt else None}


def read_fasta(path):
    """{ASV: sequence}, in the file's order."""
    seqs, name = {}, None
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                m = ASV_IN_TITLE.match(line[1:])
                if not m or m[1] in seqs:
                    raise SystemExit(f"FASTA title '{line[1:]}' does not name a single new ASV.")
                name = m[1]
                seqs[name] = ''
            elif line:
                seqs[name] += line
    return seqs


def blast_columns(t):
    ident = f"{t['identity_pct']:.2f}%" if pd.notna(t['identity_pct']) else '-'
    cover = f"{t['query_cover_pct']:.0f}%" if pd.notna(t['query_cover_pct']) else '-'
    return f"{t['closest_relatives']:<45} {cover:>5} {ident:>8}  {t['accessions']}"


def shap_rank_text(t):
    if pd.isna(t['rank']):
        return 'not kept'
    return f"{t['rank']:g}" if t['mean_abs_shap'] > 0 else 'none'


def main(folder, blast_path):
    fastas = sorted(glob.glob(os.path.join(folder, 'Table_S2_ASVs_????-??-??.fasta')))
    if len(fastas) != 1:
        raise SystemExit(f"Expected one Table_S2_ASVs_<date>.fasta in {folder}, found {len(fastas)}.")
    fasta = read_fasta(fastas[0])
    asvs = pd.read_csv(os.path.join(folder, '12_table_s2_asvs.csv'))
    if list(asvs['ASV_ID']) != list(fasta) or list(asvs['sequence']) != list(fasta.values()):
        raise SystemExit(f"{fastas[0]} does not match 12_table_s2_asvs.csv. Use the two files from the same run of 12.")
    asvs = asvs.set_index('ASV_ID')
    date = re.search(r'(\d{4}-\d{2}-\d{2})\.fasta$', fastas[0])[1]
    out_path = os.path.join(folder, f'TableS2_{date}.csv')
    if os.path.exists(out_path):
        raise SystemExit(f"{out_path} already exists.")

    header, queries = parse_blast_text(blast_path)
    if header['Program'].split()[0:1] != ['BLASTN'] or '16S_ribosomal_RNA' not in header['Database']:
        raise SystemExit(f"Expected a BLASTN search of the 16S ribosomal RNA database; this file has "
                         f"program '{header['Program']}' and database '{header['Database']}'.")
    by_asv = {}
    for q in queries:
        m = ASV_IN_TITLE.match(q['title'])
        if not m or m[1] in by_asv:
            raise SystemExit(f"Query title '{q['title']}' does not name a single new ASV.")
        by_asv[m[1]] = q
    if set(by_asv) != set(fasta):
        raise SystemExit(f"The BLAST queries do not match the FASTA. Only in the FASTA: "
                         f"{sorted(set(fasta) - set(by_asv))}; only in the BLAST file: {sorted(set(by_asv) - set(fasta))}.")
    wrong_length = [a for a in fasta if by_asv[a]['length'] != len(fasta[a])]
    if wrong_length:
        raise SystemExit(f"Query lengths differ from the FASTA for {wrong_length}.")

    rows = []
    for asv in fasta:  # 12's order: the top 20 by rank, then the rest of the Caproiciproducens ASVs
        r = asvs.loc[asv]
        rows.append({'rank': r['rank'], 'ASV_ID': asv, 'included_as': r['included_as'],
                     'silva_taxon': r['silva_taxon'], 'genus': r['genus'], 'mean_abs_shap': r['mean_abs_shap'],
                     'n_samples': r['n_samples'], 'n_studies': r['n_studies'],
                     'pct_of_genus_reads': r['pct_of_genus_reads'],
                     'cum_pct_of_genus_reads': r['cum_pct_of_genus_reads'],
                     **summarise(by_asv[asv]['hits']),
                     'n_hits': len(by_asv[asv]['hits']), 'RID': header['RID']})
    table = pd.DataFrame(rows)
    table.to_csv(out_path, index=False)

    top = table[table['included_as'].str.startswith('SHAP top')]
    genus_set = table[table['cum_pct_of_genus_reads'].notna()].sort_values('cum_pct_of_genus_reads')
    print(f"BLAST file: {blast_path} (RID {header['RID']}; {header['Database']})")
    print(f"FASTA: {fastas[0]} ({len(fasta)} ASVs, all matched)")
    print("A hit of 99% or more is the closest known relative, not a species identification.")
    print("studies: in how many studies the ASV is found (count above 0 in any modelled sample).")

    print(f"\n--- TABLE S2, PART 1: THE TOP {len(top)} ASVs BY MEAN |SHAP| ---")
    print("-" * 132)
    for _, t in top.iterrows():
        print(f"{t['rank']:>4g}  {t['ASV_ID']:<9} {t['silva_taxon']:<38} {t['n_studies']:>2}  {blast_columns(t)}")
    print("-" * 132)

    if len(genus_set):
        genus = genus_set['genus'].iloc[0]
        print(f"\n--- TABLE S2, PART 2: {genus.upper()}, THE {len(genus_set)} ASVs THAT TOGETHER HOLD "
              f"{genus_set['cum_pct_of_genus_reads'].iloc[-1]:.1f}% OF ITS READS ---")
        print(f"{'ASV':<9} {'% of genus':>10} {'cumulative':>10} {'SHAP rank':>10} {'studies':>8}  closest relative(s)")
        print("-" * 132)
        for _, t in genus_set.iterrows():
            print(f"{t['ASV_ID']:<9} {t['pct_of_genus_reads']:>9.1f}% {t['cum_pct_of_genus_reads']:>9.1f}% "
                  f"{shap_rank_text(t):>10} {t['n_studies']:>8}  {blast_columns(t)}")
        print("-" * 132)
        print("SHAP rank 'none': the ASV gets no attribution in the fit; 'not kept': removed by the prevalence filter.")
        print(f"\n{genus} reads in these {len(genus_set)} ASVs, by closest relative:")
        by_relative = (genus_set.groupby(['closest_relatives', 'at_least_99'])
                       .agg(n_asvs=('ASV_ID', 'size'), pct=('pct_of_genus_reads', 'sum'))
                       .sort_values('pct', ascending=False))
        for (relatives, at_least_99), g in by_relative.iterrows():
            n = int(g['n_asvs'])
            print(f"  {relatives}, identity {'99% or more' if at_least_99 else 'below 99%'}: "
                  f"{n} ASV{'s' if n != 1 else ''}, {g['pct']:.1f}% of the genus's reads")

    print(f"\nAt least 99% identity: {int(table['at_least_99'].sum())} of {len(table)} ASVs; "
          f"ties listed for {int((table['n_species_tied'] > 1).sum())}.")
    print(f"  Saved: {out_path}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python 12b_table_s2_from_blast.py <folder with 12's outputs> <BLAST text file>")
    main(sys.argv[1], sys.argv[2])
