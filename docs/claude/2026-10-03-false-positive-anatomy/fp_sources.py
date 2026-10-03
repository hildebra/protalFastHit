import pandas as pd, numpy as np, sys, warnings
warnings.filterwarnings('ignore')
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40); pd.set_option('display.max_rows', 300)
base = '/mnt/c/Users/hildebra/Documents/locDev/protal/local/v5'
key = ['meta_sample', 'taxon']
lin = ['meta_lineage_genus', 'meta_lineage_family', 'meta_lineage_order', 'meta_lineage_class', 'meta_lineage_phylum']
tr = pd.read_csv(f'{base}/training/training_data.tsv', sep='\t', low_memory=False)
pr = pd.read_csv(f'{base}/trained_model.predictions.tsv.gz', sep='\t', low_memory=False)
tr = tr.merge(pr[key + ['p_species'] + lin], on=key)
te = pd.read_csv(f'{base}/test/training_data.tsv', sep='\t', low_memory=False)
pte = pd.read_csv(f'{base}/trained_model.test_predictions.tsv.gz', sep='\t', low_memory=False)
te = te.merge(pte[key + ['p'] + lin], on=key).rename(columns={'p': 'p_species'})
df = pd.concat([tr.assign(set='train'), te.assign(set='test')]).reset_index(drop=True)
df['call'] = df.p_species >= 0.5
df['present'] = df.truth == 1
df['FP'] = ~df.present & df.call
df['nov'] = df.meta_novel_level.fillna('db')
df['ani_est'] = 1 - df.excess_median
one = df[df.fragments == 1].copy()
one['grp'] = np.where(one.FP, 'FP', np.where(one.present & one.call, 'TP', np.where(one.present, 'FN', 'TN')))
print('-- 1-fragment: mean_mapq bins x grp'); print(pd.crosstab(pd.cut(one.mean_mapq, [-1, 3, 10, 20, 40, 60, 80, 100, 200]), one.grp))
print('\n-- 1-fragment: identity bins x grp'); print(pd.crosstab(pd.cut(one.identity, [0, 0.93, 0.95, 0.97, 0.98, 0.99, 0.995, 1.001]), one.grp))
print('\n-- 1-fragment present by meta_rep_genome x grp'); print(pd.crosstab(one[one.present].meta_rep_genome, one[one.present].grp))
print('\n-- 1-fragment rows with identity>=0.99 by depth x grp'); hi = one[one.identity >= 0.99]; print(pd.crosstab(hi.meta_read_pairs, hi.grp))
print('\n-- 1-fragment FP: relative rank x nov'); print(pd.crosstab(one[one.FP].meta_relative_rank, one[one.FP].nov, margins=True))

# source candidates for FP whose closest simulated species is in the DB but another genus
print('\n==== FP from a present DB species of ANOTHER genus: candidate sources = present species of the sample sharing the deepest shared rank')
rankcol = {'family': 'meta_lineage_family', 'order': 'meta_lineage_order', 'class': 'meta_lineage_class', 'phylum': 'meta_lineage_phylum'}
pres = df[df.present]
rows = []
for _, r in df[df.FP & (df.nov == 'db') & df.meta_relative_rank.isin(rankcol)].iterrows():
    col = rankcol[r.meta_relative_rank]
    cand = pres[(pres.meta_sample == r.meta_sample) & (pres[col] == r[col]) & (pres.meta_lineage_genus != r.meta_lineage_genus)]
    for _, c in cand.iterrows():
        rows.append(dict(fp=r.taxon_name, fp_frag=r.fragments, fp_identity=round(r.identity, 3), rank=r.meta_relative_rank, source=c.taxon_name, source_frag=c.fragments, source_rep=c.meta_rep_genome, sample=r.meta_sample, n_cand=len(cand)))
src = pd.DataFrame(rows)
print('FP rows with candidates:', src[['fp', 'sample']].drop_duplicates().shape[0], '; of these with exactly one candidate:', (src.n_cand == 1).sum())
pairs = src[src.n_cand == 1].groupby(['fp', 'source']).agg(n=('sample', 'nunique'), fp_identity=('fp_identity', 'median'), fp_frag=('fp_frag', 'median'), source_frag=('source_frag', 'median'), rank=('rank', 'first')).sort_values('n', ascending=False)
print('\n-- (FP taxon, single candidate source) pairs (top 40)'); print(pairs.head(40))
print('\n-- how often the same source species appears for a given FP taxon across its FP samples')
multi = src.groupby('fp').agg(n_samples=('sample', 'nunique'), n_sources=('source', 'nunique'), top_source_share=('source', lambda s: s.value_counts().iloc[0] / len(s)))
print(multi[multi.n_samples >= 2].sort_values('n_samples', ascending=False).head(30))
print('\n-- source_frag distribution for single-candidate pairs (is the source abundant?)'); print(src[src.n_cand == 1].source_frag.describe())
print('\n-- for the top pairs: samples with the source present, in how many does the FP taxon have records, and is it called?')
absent = df[~df.present]
out = []
for (fpn, so), row in pairs.head(25).iterrows():
    s_pres = set(pres[pres.taxon_name == so].meta_sample)
    fp_abs = absent[(absent.taxon_name == fpn) & absent.meta_sample.isin(s_pres)]
    out.append(dict(fp=fpn, source=so, samples_source_present=len(s_pres), fp_has_records=len(fp_abs), fp_called=int(fp_abs.FP.sum()), med_identity=round(fp_abs.identity.median(), 3) if len(fp_abs) else None, med_frag=fp_abs.fragments.median() if len(fp_abs) else None))
print(pd.DataFrame(out))

# gene-pattern separability for >=3 fragments
print('\n==== rows with >=3 fragments: FP vs TP vs FN vs TN(identity>=0.95) feature quartiles')
m = df[df.fragments >= 3].copy()
m['grp'] = np.where(m.FP, 'FP', np.where(m.present & m.call, 'TP', np.where(m.present, 'FN', 'TN')))
m = m[(m.grp != 'TN') | (m.identity >= 0.95)]
cols = ['fragments', 'identity', 'top_identity', 'ani_est', 'excess_high_share', 'low_identity_share', 'hit_gene_fraction', 'gene_presence_ratio', 'depth_cv', 'gene_dispersion', 'conserved_hit_share', 'conserved_fast_depth_ratio', 'conserved_fast_record_ratio', 'conserved_fast_kept_ratio', 'mean_mapq', 'low_mapq_share', 'congener_fit_share', 'genus_skew', 'relative_distance', 'relative_close_share', 'em_own_share', 'variant_sites_per_kb', 'p_species']
for c in cols:
    q = m.groupby('grp')[c].quantile([0.25, 0.5, 0.75]).unstack()
    print(f'{c:30s}', '  '.join(f'{k}: {q.loc[k,0.25]:.3g}/{q.loc[k,0.5]:.3g}/{q.loc[k,0.75]:.3g}' for k in ['TP', 'FN', 'FP', 'TN']))
print(m.grp.value_counts())
print('\n-- >=3 fragments FP: relative rank x nov'); print(pd.crosstab(m[m.FP].meta_relative_rank, m[m.FP].nov, margins=True))
print('\n-- >=3 fragments FP by identity bin x nov'); print(pd.crosstab(pd.cut(m[m.FP].identity, [0, 0.95, 0.96, 0.97, 0.98, 0.99, 1.001]), m[m.FP].nov, margins=True))
print('-- >=3 fragments TP by identity bin'); print(pd.cut(m[m.grp == 'TP'].identity, [0, 0.95, 0.96, 0.97, 0.98, 0.99, 1.001]).value_counts().sort_index())
print('-- >=3 fragments FN by identity bin'); print(pd.cut(m[m.grp == 'FN'].identity, [0, 0.95, 0.96, 0.97, 0.98, 0.99, 1.001]).value_counts().sort_index())
