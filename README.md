# Natural statistics of host odours predict species-specific olfactory behaviours in drosophilids

Analysis code and data for Gong *et al.* (Current Biology, in revision).

The paper ranks the volatile compounds of each species' host fruits by how much they vary across fruits (the **PC score**), and tests whether high-ranking compounds are behaviourally relevant for *Drosophila melanogaster* (generalist), *D. sechellia* (noni specialist) and *D. erecta* (pandan specialist). This repository contains the code for the PC-score analyses, their robustness checks and the host-versus-non-host random forests.

## Contents

| Folder | What it does | Figure panels |
|---|---|---|
| `pc_on_averages/For_main_figures` | Within-fruit z-scoring, PCA and PC score for the three datasets; PC score by behavioural class for *D. melanogaster* | Fig. 1b–d, 1f, 2a; Fig. S1a |
| `pc_on_averages/CLR` | Centred log-ratio (CLR) normalisation as an alternative to within-fruit z-scoring (*D. melanogaster*) | Fig. S1b |
| `noni_ripeness_projection` | PCA of laboratory noni fruits staged for ripeness, and projection of the field-collected fruits onto it; used to assign field ripeness categories | Fig. 1j |
| `resilience` | Resilience of PC-score rankings to the fruits sampled: leave-one-fruit-out (noni, pandan), leave-three-fruits-out (*D. melanogaster*) and a citrus-exclusion control | Fig. S1d, S1e, S1h |
| `robustness` | Robustness of PC scores to replicate selection (one random replicate per fruit, 1,000 combinations) | Fig. 1g, 2b |
| `random_forest` | Host-versus-non-host random forests and their comparison with the PC-score rankings | Fig. 1m–p, 2h–k |

The ripeness and fruit-identity classifiers (Fig. 1k–l, 2f–g, S1f–g and S1i–l) are in a separate repository; see [Related code](#related-code).

## Requirements

- **MATLAB** R2020a or later, with the Statistics and Machine Learning Toolbox (`pca`, `boxplot`, `tiedrank`). `rf_vs_pca_full.m` needs base MATLAB only. The Fig. 1j script also needs the Curve Fitting Toolbox (`cscvn`, `fnplt`) and the `rgb` function from the MATLAB File Exchange.
- **Python** 3.9 or later, with `pandas`, `numpy`, `scikit-learn` and `openpyxl`. The published results used scikit-learn 1.7.1; other versions can give slightly different importance values.

## How to run

Each folder is self-contained: its input files sit next to the script. Run each script from within its own folder.

### `pc_on_averages/For_main_figures`

`chem_pca_plus_PCscore_clean_v3.m` computes the PC score of every compound. Select the dataset by editing the `inputFile` line:

- `dmel_data.csv`: *D. melanogaster* host fruits. Only for this dataset, the script also reads `behaviour-anotations.xlsx` and plots PC score by behavioural class (Fig. 1d).
- `dsec_data_fruit_avgs.csv`: noni fruits.
- `pandanfruit_averages.csv`: pandan fruits.

The script also plots cumulative variance explained, variance-weighted loadings and the fruits in PC space, and prints the most positive and negative loadings on PC1–PC3.

### `pc_on_averages/CLR`

1. `transform_CLR/make_CLR.m` reads `dmel_data.csv` and applies the CLR transform. Before the transform, negative values are set to zero and zeros are replaced with half the smallest positive value in that fruit. It writes `Dmel_CLR_transformed.csv` (fruits × compounds, with compounds labelled `Odour_1`, `Odour_2`, …) and `Dmel_odour_name_key.csv` (the compound name for each label).
2. `run_CLR/PCA_plus_plot.m` runs the PCA on the CLR data, computes the PC score, plots it by behavioural class and writes `Annotated_CLR_PC_scores_from_MATLAB.csv`. It needs the two files from step 1 and `behaviour-anotations.xlsx` in its folder.

### `noni_ripeness_projection`

`noniRipenessPCAprojectioncolors.m` places the field-collected noni fruits on the ripeness axis of laboratory fruits staged for ripeness (Fig. 1j). It uses the 14 odorants shared by the two datasets and reads them from the `input_for_matlab` sheet of `input.xlsx`. The script:

1. converts each laboratory sample to relative abundances, z-scores it within the sample, and runs a PCA on the laboratory samples;
2. normalises each field sample the same way, standardises it with the laboratory mean and standard deviation of each odorant, and projects it onto the laboratory principal components;
3. plots the laboratory and field samples in PC space, with a curve through the mean of each laboratory ripeness stage, and plots each fruit's mean position along PC1;
4. plots the variance-weighted PC loadings and prints the odorants with the largest loadings on PC1–PC3.

PC1 orders the laboratory fruits by ripeness (green → brown), and its strongest loadings are methyl hexanoate (less ripe) and hexanoic and octanoic acid (more ripe). Field fruits were assigned to four ripeness categories from their positions along this axis. One laboratory sample (G2) contains none of the 14 shared odorants, so it is undefined after normalisation and is excluded from the PCA (14 laboratory samples are used). Field fruit S11 is labelled "Field 11a" in the script.

### `resilience`

- `pcscore_leave_fruit_out.m` recomputes the PC score leaving out each fruit in turn, for noni (12 iterations) and pandan (8 iterations).
- `pcscore_leave3_dmel.m` does the same for *D. melanogaster*, removing three fruits at a time. It uses 1,000 unique random combinations (fixed seed) and displays the 10 top-ranked compounds.
- `pcscore_dmel_no_citrus.m` removes the three citrus fruits (grapefruit, mandarin and orange) and reports how the rank and score of D-limonene change. It then repeats the leave-three-out analysis on the remaining 30 fruits.

### `robustness`

`pc_score_robustness.m` draws 1,000 random combinations of one replicate per fruit (fixed seed; duplicate combinations removed), computes the PC score for each and plots the distributions. Set the input file, the output file names and `fruit_name` at the top of the script for each dataset:

- noni: `input_Dsec_robustness.csv`
- pandan: `input_Dere_robustness.csv`

### `random_forest`

1. `random_forest_figures_used/rf_host_vs_nonhost.py` trains random forests to separate noni (or pandan) samples from all other fruits:

   ```bash
   python rf_host_vs_nonhost.py
   ```

   It uses two feature sets:

   - **Permissive** (75 compounds): all compounds detected in noni or pandan.
   - **Conservative** (44 compounds): only those also present in the *D. melanogaster* dataset, so that a compound absent there for methodological reasons cannot drive the classification.

   Compounds are matched across datasets by InChIKey and z-scored within each sample. Fruits are held out one at a time (leave-one-group-out cross-validation), and importances are averaged across folds and 25 random seeds. The script writes four ranking files (`{noni,pandan}_vs_others_RF_named_{conservative,permissive}.xlsx`) and prints the accuracy and top six compounds for each.

2. `Random_forest_PCA_comparison/rf_vs_pca_full.m` takes the four ranking files and the PC-score rankings (`noni-pcscpre-inkey_named.xlsx`, `pandan-pcscore-inkey_named.xlsx`). It keeps the compounds common to both, compares the top six from each, and draws the alluvial (Sankey) diagrams and hypergeometric overlap tests as PDF files. The top-N size is set by `N` at the top of the script.

## Data

| File | Contents |
|---|---|
| `dmel_data.csv` | Published *D. melanogaster* host-fruit volatiles (Dweck *et al.*, 2018, *Cell Reports*): 2,491 compounds × 33 fruits, relative abundances |
| `dsec_data_fruit_avgs.csv` | Field-collected noni fruits (this study): 38 compounds × 12 fruits; mean of 4–6 replicate samples per fruit |
| `pandanfruit_averages.csv` | Field-collected pandan fruits (this study): 45 compounds × 8 fruits; mean of 3–4 replicate samples per fruit |
| `input_Dsec_robustness.csv`, `input_Dere_robustness.csv` | Replicate-level noni (55 samples) and pandan (28 samples) data |
| `noni_fruits.xlsx`, `Pandan_fruits.xlsx`, `dmel_data__1_.xlsx` | The same data with CAS numbers, InChI and InChIKeys, used to match compounds across datasets |
| `inchi_odour_mapping.csv` | Compound names for InChIKeys |
| `noni-pcscpre-inkey_named.xlsx`, `pandan-pcscore-inkey_named.xlsx` | PC scores (all fruits) with InChIKeys |
| `input.xlsx` | Noni ripeness data. Sheet `input_for_matlab`: the 14 odorants shared by the laboratory and field datasets, for 55 field replicate samples (columns 2–56) and 15 laboratory samples (columns 57–71: green G, green-yellow GY, yellow Y, yellow-brown YB and brown B, three fruits each). Sheet `field_renamed`: field fruit averages. Sheet `auer_raw`: the full laboratory dataset (Auer *et al.*, 2020, *Nature*), 46 compounds × 15 fruits |
| `behaviour-anotations.xlsx` | Literature behavioural annotations for *D. melanogaster*: 1 attractive, 0 neutral, 0.5 unclear (conflicting reports), −1 aversive (excluded from plots); blank = not annotated |

For the field data, integrated GC–MS peak areas below 10,000 were set to zero before averaging. In all data files, compounds are rows and fruits (or samples) are columns, with compound names in the first column.

## PC score

Each fruit's profile is z-scored across its compounds (within-fruit z-scoring). PCA is then run with fruits as observations and compounds as variables. The PC score of compound *j* is

PC score*ⱼ* = Σ<sub>*i*=1..3</sub> |loading<sub>*ij*</sub>| × (% variance explained by PC*ᵢ*)

## Related code

- **Ripeness and fruit-identity classifiers** (noni ripeness regression, pandan fruit-identity classification and ripeness prediction on independently staged fruit): [stootoon/drosophila-natural-odour-stats-release](https://github.com/stootoon/drosophila-natural-odour-stats-release), notebooks `sechellia-demo.ipynb` (Fig. 1k–l, S1f–g), `erecta-demo.ipynb` (Fig. 2f–g, S1i–j) and `ripeness-demo.ipynb` (Fig. S1k–l). That repository also contains `noni-random-forest.ipynb`, the original single-seed implementation of the conservative random forest that `rf_host_vs_nonhost.py` extends.
- **Oviposition assays**: [PrietoGodinoLab/High-throughput_drosophila_oviposition_preference_assay](https://github.com/PrietoGodinoLab/High-throughput_drosophila_oviposition_preference_assay)

## Citation and contact

If you use this code or data, please cite Gong *et al.* (Current Biology, in revision; full reference to be added on publication).

Contact: Lucia L. Prieto-Godino (lucia.prietogodino@crick.ac.uk)
