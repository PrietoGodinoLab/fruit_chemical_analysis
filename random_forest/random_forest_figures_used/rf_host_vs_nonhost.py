"""
rf_host_vs_nonhost.py
Random-forest host-vs-non-host odour analysis (noni / pandan), with two feature sets:

  CONSERVATIVE (44 odours): (noni InChIKeys UNION pandan InChIKeys) INTERSECT dmel InChIKeys.
      Reproduces the collaborator's pipeline. Excludes odours (e.g. free acids) not
      detected in the SPME dmel dataset, where a non-detection is ambiguous.

  PERMISSIVE   (75 odours): noni InChIKeys UNION pandan InChIKeys (all PDMS-detected
      odours), regardless of the dmel dataset; dmel non-detections are treated as
      true absences (value 0). Adds back the free acids.

Pipeline:
  - pool molecules by InChIKey; drop dmel 'unsure' InChIKeys; sum duplicate InChIKeys
  - combine all samples on the chosen feature set; missing molecule -> 0
  - z-score per sample (per column/observation) across the feature set
  - group samples by fruit; label by task (noni-vs-non-noni or pandan-vs-non-pandan)
  - RandomForestClassifier (scikit-learn defaults) under LeaveOneGroupOut CV
  - one accuracy + one importance vector per held-out fold; average importances across folds
Importances here are additionally averaged over several random seeds for stability.

Requires: pandas, numpy, scikit-learn, openpyxl.
"""
import pandas as pd, numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneGroupOut

META = ['Odour','CAS','InChI','InChIKey']
FILES = {'dmel':'dmel_data__1_.xlsx','noni':'noni_fruits.xlsx','pandan':'Pandan_fruits.xlsx'}

def load_long(f, dataset):
    df = pd.read_excel(f)
    df['InChIKey'] = df['InChIKey'].astype(str).str.strip()
    df = df[~df['InChIKey'].isin(['unsure','nan',''])]
    samp = [c for c in df.columns if c not in META]
    long = df.melt(id_vars=['InChIKey'], value_vars=samp, var_name='sample', value_name='val')
    long['val'] = pd.to_numeric(long['val'], errors='coerce').fillna(0.0)
    long = long.groupby(['InChIKey','sample'], as_index=False)['val'].sum()   # sum duplicate InChIKeys
    if dataset == 'dmel':
        long['group'] = 'mel_' + long['sample']
    else:
        long['group'] = dataset + '_' + long['sample'].str.split('_rep').str[0]
    long['sample'] = dataset + '::' + long['sample']
    long['dataset'] = dataset
    return long

def build(alld, feature_keys):
    sub = alld[alld['InChIKey'].isin(feature_keys)]
    X = sub.pivot_table(index='sample', columns='InChIKey', values='val', aggfunc='sum', fill_value=0.0)
    X = X.reindex(columns=feature_keys, fill_value=0.0)
    meta = alld[['sample','group','dataset']].drop_duplicates().set_index('sample').reindex(X.index)
    mu = X.mean(axis=1).values[:,None]
    sd = X.std(axis=1, ddof=0).values[:,None].copy(); sd[sd==0] = 1.0
    Xz = (X.values - mu) / sd
    return Xz, meta['group'].to_numpy(), meta['dataset'].to_numpy(), list(X.columns)

def run_rf(alld, feature_keys, focal, seed):
    Xz, groups, datasets, cols = build(alld, feature_keys)
    y = (datasets == focal).astype(int)
    logo = LeaveOneGroupOut(); imps=[]; accs=[]
    for tr, te in logo.split(Xz, y, groups):
        clf = RandomForestClassifier(random_state=seed)
        clf.fit(Xz[tr], y[tr])
        accs.append((clf.predict(Xz[te]) == y[te]).mean())
        imps.append(clf.feature_importances_)
    return pd.Series(np.mean(imps,0), index=cols), float(np.mean(accs))

def avg_run(alld, feature_keys, focal, seeds=range(25)):
    mis=[]; accs=[]
    for s in seeds:
        mi, acc = run_rf(alld, feature_keys, focal, s); mis.append(mi); accs.append(acc)
    return pd.concat(mis, axis=1).mean(axis=1), float(np.mean(accs))

def name_map():
    nm={}
    for f in FILES.values():
        df = pd.read_excel(f, usecols=['Odour','InChIKey'])
        for _,r in df.iterrows():
            k = str(r['InChIKey']).strip()
            if k not in ('unsure','nan','') and k not in nm and isinstance(r['Odour'], str):
                nm[k] = r['Odour'].replace('\xa0',' ').strip()
    try:
        mp = pd.read_csv('inchi_odour_mapping.csv')
        for _,r in mp.iterrows(): nm[str(r['InChiKey']).strip()] = str(r['Odour']).strip()
    except FileNotFoundError:
        pass
    return nm

if __name__ == '__main__':
    alld = pd.concat([load_long(FILES[d], d) for d in FILES], ignore_index=True)
    nk = set(alld[alld.dataset=='noni'].InChIKey); pk = set(alld[alld.dataset=='pandan'].InChIKey)
    dk = set(alld[alld.dataset=='dmel'].InChIKey)
    CONS = sorted((nk | pk) & dk)        # 44
    PERM = sorted(nk | pk)               # 75
    nm = name_map(); nn = lambda k: nm.get(k, k)
    SEEDS = range(25)                    # raise for more stability; 8 is enough for the top odours

    for focal in ['noni','pandan']:
        for tag, feats in [('conservative',CONS), ('permissive',PERM)]:
            imp, acc = avg_run(alld, feats, focal, SEEDS)
            df = pd.DataFrame({'InChIKey':imp.index,
                               'chemical_name':[nn(k) for k in imp.index],
                               'mean':imp.values}).sort_values('mean', ascending=False)
            df.to_excel(f'{focal}_vs_others_RF_named_{tag}.xlsx', index=False)
            print(f'{focal} {tag}: accuracy={acc:.3f}; top6=' +
                  ', '.join(df.head(6)['chemical_name']))
