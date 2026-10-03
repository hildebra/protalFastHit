import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
pd.set_option('display.width', 250)
base = '/mnt/c/Users/hildebra/Documents/locDev/protal/local/v5'
key = ['meta_sample', 'taxon']
tr = pd.read_csv(f'{base}/training/training_data.tsv', sep='\t', low_memory=False)
pr = pd.read_csv(f'{base}/trained_model.predictions.tsv.gz', sep='\t', low_memory=False)
tr = tr.merge(pr[key + ['p_species']], on=key)
te = pd.read_csv(f'{base}/test/training_data.tsv', sep='\t', low_memory=False)
pte = pd.read_csv(f'{base}/trained_model.test_predictions.tsv.gz', sep='\t', low_memory=False)
te = te.merge(pte[key + ['p']], on=key).rename(columns={'p': 'p_species'})
df = pd.concat([tr, te]).reset_index(drop=True)
df['present'] = df.truth == 1
df['call'] = df.p_species >= 0.5
df['grp'] = np.where(~df.present & df.call, 'FP', np.where(df.present & df.call, 'TP', np.where(df.present, 'FN', 'TN')))
for n, sel in [('1 fragment', df.fragments == 1), ('2 fragments', df.fragments == 2), ('3-9', (df.fragments >= 3) & (df.fragments < 10)), ('>=10', df.fragments >= 10)]:
    d = df[sel]
    print(f'\n== {n}: linked_share bins x group')
    print(pd.crosstab(pd.cut(d.linked_share, [-0.01, 0.0, 0.25, 0.5, 0.75, 0.999, 1.0]), d.grp))
one = df[df.fragments == 1]
print('\n== 1 fragment, identity>=0.97: linked_share==0 vs ==1 by group')
d = one[one.identity >= 0.97]
print(pd.crosstab(d.linked_share.round(2), d.grp))
print('\n== 1 fragment, by read length (meta_read_length): share with linked_share==1 among TP and FP')
for rl, g in one.groupby('meta_read_length'):
    print(rl, {k: f"{(g[g.grp==k].linked_share==1).mean():.2f} (n={int((g.grp==k).sum())})" for k in ['TP', 'FP', 'FN', 'TN']})
print('\n== se table: linked_share values (should be 0)')
se = pd.read_csv(f'{base}/training/training_data_se.tsv', sep='\t', usecols=['linked_share', 'fragments'], low_memory=False)
print(se.linked_share.describe())
