import pandas as pd, numpy as np, sys
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40); pd.set_option('display.max_rows', 200)
base = sys.argv[1] if len(sys.argv) > 1 else '/mnt/c/Users/hildebra/Documents/locDev/protal/local/v5'
rt = sys.argv[2] if len(sys.argv) > 2 else ''   # '', '_se', '_pb', '_ont'
sfx = rt
key = ['meta_sample', 'taxon']
tr = pd.read_csv(f'{base}/training/training_data{sfx}.tsv', sep='\t', low_memory=False)
pr = pd.read_csv(f'{base}/trained_model{sfx}.predictions.tsv.gz', sep='\t', low_memory=False)
tr = tr.merge(pr[key + ['p_species', 'meta_lineage_genus']], on=key)
te = pd.read_csv(f'{base}/test/training_data{sfx}.tsv', sep='\t', low_memory=False)
pte = pd.read_csv(f'{base}/trained_model{sfx}.test_predictions.tsv.gz', sep='\t', low_memory=False)
te = te.merge(pte[key + ['p', 'meta_lineage_genus']], on=key).rename(columns={'p': 'p_species'})

def fbin(f):
    return pd.cut(f, [0, 1, 2, 9, 99, 999, 1e12], labels=['1', '2', '3-9', '10-99', '100-999', '>=1000'])
def ibin(i):
    return pd.cut(i, [0, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0001], labels=['<0.93', '0.93-0.95', '0.95-0.96', '0.96-0.97', '0.97-0.98', '0.98-0.99', '>=0.99'])

for name, df in [('TRAINING (p with species held out)', tr), ('TEST', te)]:
    df = df.copy()
    df['call'] = df.p_species >= 0.5
    df['FP'] = (df.truth == 0) & df.call
    df['FN'] = (df.truth == 1) & ~df.call
    df['present'] = df.truth == 1
    df['ani_est'] = 1 - df.excess_median
    df['nov'] = df.meta_novel_level.fillna('db-species-closest')
    print(f'\n===== {name}: rows {len(df)}, present {df.present.sum()}, FP {df.FP.sum()}, FN {df.FN.sum()}')
    print('\n-- FP by the deepest rank shared with any simulated species (meta_relative_rank) x whether that closest species is held out (meta_novel_level)')
    fp = df[df.FP]
    print(pd.crosstab(fp.meta_relative_rank, fp.nov, margins=True))
    print('\n-- all ABSENT rows by the same, for the FP rate')
    ab = df[~df.present]
    print(pd.crosstab(ab.meta_relative_rank, ab.nov, margins=True))
    print('\n-- FP by fragments bin x identity bin (median read identity)')
    print(pd.crosstab(fbin(fp.fragments), ibin(fp.identity), margins=True))
    print('\n-- PRESENT (all) by fragments bin x identity bin')
    pres = df[df.present]
    print(pd.crosstab(fbin(pres.fragments), ibin(pres.identity), margins=True))
    print('\n-- PRESENT called (TP) by fragments x identity')
    print(pd.crosstab(fbin(pres[pres.call].fragments), ibin(pres[pres.call].identity), margins=True))
    print('\n-- FP by fragments bin x ani_est bin (1 - excess_median)')
    print(pd.crosstab(fbin(fp.fragments), ibin(fp.ani_est), margins=True))
    print('\n-- PRESENT by fragments bin x ani_est bin')
    print(pd.crosstab(fbin(pres.fragments), ibin(pres.ani_est), margins=True))
    print('\n-- medians per class (absent split by source)')
    df['cls'] = np.where(df.present, 'present', np.where(df.nov == 'db-species-closest', 'absent, closest species in db', 'absent, closest species held out'))
    cols = ['fragments', 'identity', 'top_identity', 'ani_est', 'excess_high_share', 'low_identity_share', 'hit_gene_fraction', 'conserved_hit_share', 'conserved_fast_record_ratio', 'genus_skew', 'relative_distance', 'em_own_share', 'congener_fit_share', 'mean_mapq', 'low_mapq_share']
    print(df.groupby('cls')[cols].median().T)
    print('\n-- the same, FP only vs TP only, fragments < 10')
    small = df[df.fragments < 10]
    g = small.assign(grp=np.where(small.FP, 'FP', np.where(small.present & small.call, 'TP', np.where(small.FN, 'FN', 'TN')))).groupby('grp')
    print(g[cols].median().T)
    print(g.size())
    print('\n-- rule grid: veto rows with identity < t (and ani_est < t2) among fragments < N: FP removed / TP removed')
    for N in [3, 10, 1e9]:
        for t in [0.93, 0.94, 0.95, 0.96, 0.97]:
            m = (df.fragments < N) & (df.identity < t)
            print(f'  N<{N:g} identity<{t}: FP removed {int((m & df.FP).sum()):4d} of {int(df.FP.sum())}, TP removed {int((m & df.present & df.call).sum()):4d} of {int((df.present & df.call).sum())}')
        for t in [0.93, 0.94, 0.95, 0.96]:
            m = (df.fragments < N) & (df.ani_est < t)
            print(f'  N<{N:g} ani_est<{t}: FP removed {int((m & df.FP).sum()):4d}, TP removed {int((m & df.present & df.call).sum()):4d}')
    # genus-level spray pattern
    print('\n-- genus pattern: per (sample, genus): species with records, top share; for FP rows vs TP rows vs TN rows with records')
    df['genus'] = df.meta_lineage_genus
    gg = df.groupby(['meta_sample', 'genus'])
    df['genus_n_species_with_records'] = gg.fragments.transform('size')
    df['genus_fragments'] = gg.fragments.transform('sum')
    df['genus_top_share'] = gg.fragments.transform('max') / df.genus_fragments
    df['own_share'] = df.fragments / df.genus_fragments
    df['genus_present_species'] = gg.present.transform('sum')
    df['grp'] = np.where(df.FP, 'FP', np.where(df.present & df.call, 'TP', np.where(df.FN, 'FN', 'TN')))
    print(df.groupby('grp')[['genus_n_species_with_records', 'genus_top_share', 'own_share', 'genus_present_species']].describe().T.round(3))
    print('\n-- FP: is the FP the top species of its genus in the sample? and does its genus hold a present species?')
    fp = df[df.FP]
    print(pd.crosstab(fp.own_share == fp.genus_top_share, fp.genus_present_species > 0, margins=True))
    print('\n-- FP with no present species in the genus (a held-out species or distant reads): siblings with records')
    print(fp[fp.genus_present_species == 0].genus_n_species_with_records.describe())
    print('\n-- TP with fragments<10: siblings with records')
    print(df[(df.grp == 'TP') & (df.fragments < 10)].genus_n_species_with_records.describe())
    print('\n-- FP by depth (meta_read_pairs)')
    print(df[df.FP].groupby('meta_read_pairs').size())
    print(df.groupby('meta_read_pairs').apply(lambda d: pd.Series({'absent': int((~d.present).sum()), 'absent_with_1frag': int(((~d.present) & (d.fragments == 1)).sum()), 'FP': int(d.FP.sum()), 'FP_1frag': int((d.FP & (d.fragments == 1)).sum()), 'present': int(d.present.sum()), 'FN': int(d.FN.sum())})))
