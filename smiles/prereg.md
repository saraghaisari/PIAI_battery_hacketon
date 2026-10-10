# Pre-registration: does molecular structure add anything beyond physics?

Written before any script for this study exists and before any model has been trained on these feature
sets. Authors: Sara Ghaisari (design), with Claude Code. Repository state at writing: Gelavizh Ahmadi's
history up to `d6e9f52`.

What we had already seen when writing this: the data summary below (row counts only, no targets
by feature set), Gelavizh's published fold A/B results for her own models, and our rerun of her `kan.py`
and `final.py`. Nothing in this study's feature sets has been fitted.

## Research question

Does a solvent's molecular structure, encoded from SMILES, improve conductivity prediction for a solvent
the model has never seen, beyond what bulk physical properties already provide?

## Hypotheses

- **H1.** Physics features (set A) beat SMILES descriptors (set B) on unseen solvents.
- **H2a.** Using both (set C) is no worse than physics alone (set A).
- **H2b.** Set C beats set A on solvents whose bulk properties are unusual or poorly known (subset U, below).
- **H0 (null).** Adding SMILES features to physics (C vs A) changes the error by less than the spread
  across random seeds.

Any outcome is reportable; failed hypotheses are reported as failed.

## Data facts that shaped the design (Step 0 of this study, counts only)

- Train: 6759 rows, 37 solvents, 87 distinct mixtures, 12 salts. Test: DEC, DMC, EC, PC only;
  salts LiBF4, LiPF6, LiBOB.
- Gelavizh's property table (`PROPS` in `pipeline.py`: M, density, ε, viscosity at about 25 °C) covers 24
  solvents and **6565 / 6759 train rows**. The 14 solvents without values each have fewer than 50 rows.
- Coverage is very uneven. Solvents with ≥ 50 covered rows (rows / distinct sources): PC 4579/12,
  EC 3506/16, EMC 1309/8, EA 710/1, DME 541/5, TFP 540/1, DMC 396/6, Sulfolane 253/1, 2-MeTHF 236/1,
  3-Glyme 171/1, Methylene chloride 164/1, 2-Glyme 136/2, MA 120/1, THF 85/1, DMSO 76/2.
- **Confound:** 9 of these 15 solvents come from a single publication. Holding out such a solvent also
  holds out a lab. Leave-one-solvent-out error therefore mixes "new solvent" with "new source". We
  cannot separate them and will state this as a limitation.

## Feature sets

All sets share the same non-solvent features: 1000/T, molality (√m, m, m²; mol/L rows converted with
`mixture_density_g_cm3` as in `pipeline.py`), a mol/L-unit flag, and salt one-hot. Mixture weighting
uses `solvent_fracs_mol`.

| Set | Solvent information |
|---|---|
| A. Physics | Mole-fraction-weighted M, ρ, ε, ln η (`PROPS` values). Donor number not used: it is not in the table, and we will not add it now. |
| B. SMILES | Mole-fraction-weighted RDKit descriptors, fixed list of 18: MolWt, MolMR, LabuteASA, TPSA, MolLogP, NumHAcceptors, NumHDonors, NumRotatableBonds, RingCount, NumAromaticRings, FractionCSP3, HeavyAtomCount, NumHeteroatoms, and atom counts of O, N, F, S, Cl. |
| B5 | 5 of the 18, chosen inside each training fold by mutual information with log_k. |
| C. Both | A + B. |
| N. Negative control | One-hot solvent identity, mole-fraction weighted (the held-out solvent's column is zero in training). |

**Property table.** `PROPS` is frozen as it stands at `d6e9f52`. Before Step 1 we check each value
against a literature source and record the source. Any correction is logged as an amendment before
training. Rows containing a solvent without properties are dropped **from every set**, so all sets
are scored on identical rows.

## Models

The same settings for every feature set; only the features change.

1. **GBM:** `HistGradientBoostingRegressor` (max_iter 500, learning_rate 0.05, max_leaf_nodes 31,
   min_samples_leaf 20, no early stopping, monotone decreasing in 1000/T). It is deterministic, so it gets one run per fold.
2. **MLP:** 2 hidden layers × 64, SiLU, Adam (lr 1e-3, weight decay 1e-4), 300 epochs, batch 256,
   inputs standardised on the training fold. 5 seeds (0 to 4).

KAN variants are exploratory only and are not part of any verdict.

## Validation

- **Primary: leave-one-solvent-out (LOSO)** over the 13 solvents with ≥ 50 covered rows **excluding
  PC and EC**. Each fold holds out every row containing that solvent. PC and EC are base solvents
  (68% and 52% of rows); holding them out removes most of the data, so their folds are reported
  separately as a stress test and are not used in any verdict.
- **Secondary:** Gelavizh's fold A (all PC+EA rows) and fold B (all EC+PC rows).
- **Final:** the Kaggle DEC test set, scored once on the leaderboard by the best set, never used for tuning.
  This session cannot reach Kaggle, so a co-author submits the file.

## Metrics

- RMSE on log_k per fold, then the **mean per-fold RMSE** over the 13 primary folds (each solvent
  weighted equally).
- Also reported: the trimmed RMSE without rows below 10⁻³ mS/cm (log_k < −3), mean bias, RMSE for T < 250 K
  and T ≥ 250 K, per salt, and per solvent.
- Seed spread σ_seed: for the MLP, the mean over primary folds of the standard deviation of fold RMSE across 5 seeds.
- Lower bound: the training-fold mean of log_k, per fold.

## Decision rules (fixed now)

Δ(X, Y) = mean per-fold RMSE of X minus that of Y. CI = 95% paired bootstrap over the 13 held-out
solvents (10,000 resamples, seed 0). Rules are applied separately for GBM and MLP (MLP uses the
seed-mean prediction), and a hypothesis is **held** only if it holds for both models.

| | Held if |
|---|---|
| H1 | Δ(B, A) > 0 and its CI lies entirely above 0 |
| H2a | Upper CI bound of Δ(C, A) < +0.02 (non-inferiority margin 0.02 log10) |
| H2b | Mean RMSE on subset U is lower for C than A by more than σ_seed (MLP); for GBM, C lower than A in ≥ 2 of the U folds |
| H0 | MLP only: \|Δ(C, A)\| < σ_seed |
| Setup check 1 | The negative control N is worse than A (Δ(N, A) > 0) |
| Setup check 2 | A beats the lower bound in ≥ 11 of 13 primary folds |

If a setup check fails, H1 to H0 are still computed but are flagged as not interpretable.

**Subset U** (fixed now from `PROPS`): primary-fold solvents whose ε or η lies outside the range of
the other primary-fold solvents plus PC and EC (computed in Step 0, before training), together with TFP,
which Gelavizh's table marks as the least certain. FEC and MOEMC are also marked uncertain but have
fewer than 50 rows, so they are not primary folds.

**Power note.** With 13 held-out solvents the bootstrap has low power, so small effects may not
reach significance. Effects are reported as estimates with CIs, not just as held/failed.

## Pitfalls checked and reported

- **Memorisation:** with 18 descriptors, each of the 24 solvents has a unique descriptor vector. LOSO
  exposes this, and the comparison with N calibrates it.
- **Unequal width:** B (18 columns) vs A (4 columns), so B5 is reported alongside B.
- **Leakage:** scaling and B5 selection are fit on training folds only.
- **Overlap between A and B:** we report Pearson correlations between every A and B column across the 24 solvents.

## Environment

Pinned in `requirements.txt`: Python 3.13, numpy 2.5.3, pandas 2.3.3, scikit-learn 1.9.1,
torch 2.14.1, scipy 1.18.1, matplotlib 3.11.2, rdkit 2026.03.6. Gelavizh's code does not run on pandas 3.

## Amendment 1 (before any training script exists; no model has been fitted)

**Reason.** Molar mass and density of every solvent can be recovered from the competition files, so
these two properties need no external table.

**Changes to set A.**
- **M:** RDKit molecular weight from `solvent_properties.csv` SMILES. This matches `PROPS` to 0.01 g/mol.
- **ρ:** recovered from `mixture_density_g_cm3` (`smiles/props_kaggle.py`). According to `metaData.csv`,
  that column is an ideal additive-volume mix of the pure densities. 183 distinct mixtures determine
  all 38 pure densities: rank 38, max residual 4.6e-9. Two values differ from `PROPS` by more than
  0.02 g/cm³: **TFP** (1.487 vs 1.59) and **MOEMC** (1.070 vs 1.10). The recovered values are the ones
  the dataset's own mole fractions were computed with, so they replace `PROPS`.
- **ε and η:** a draft for all 38 solvents is in `smiles/solvent_eta_eps_draft.csv`. These values
  are recalled from standard references (Xu 2004 Chem. Rev. 104:4303; Riddick, Bunger & Sakano 1986;
  CRC Handbook) and have **not yet been checked against the printed sources**, because this session
  cannot reach any literature site.
  - Where both exist, the draft agrees with `PROPS` within 2.5%. This is **not** an independent check:
    both likely derive from the same commonly quoted values.
  - There is no value at all for TFP, MOEMC, DMM or Propylsulfone, plus η for 3-MeSulfolane.
  - FEC, Ethyldiglyme, Pseudocumeme and TEOS are low confidence.

**Freeze rule.** ε and η are frozen only after each value has been checked against a printed or
primary source and its reference recorded. That happens in a later amendment, before Step 1. Until
then, `PROPS` ε and η remain the working values. The row-coverage rule is unchanged: a row is used
only if every solvent in it has ε and η.

**Noted for the paper.** TFP (540 rows, a single source) has no verified ε or η, and its `PROPS`
density was 7% off. Its dataset source is probably Ding, Xu & Jow 2002 (doi 10.1149/1.1513556),
which should be checked for measured values.

**Note (literature check in progress, no model fitted).** Verified values are collected in
`smiles/solvent_eta_eps_lit.csv`, one row per value with its source and table. Xu 2004's "DMM" is
dimethoxymethane (M 76), not this dataset's DMM (dipropylene glycol dimethyl ether, M 162), so it
must not be used for DMM.

## Amendment 2 (before any training script exists; no model has been fitted)

**Set A ε and η are frozen** in `smiles/solvent_props_final.csv`, built by `smiles/build_props.py`:
- A cited literature value (`smiles/solvent_eta_eps_lit.csv`) where one exists, taking the value
  measured closest to 25 °C, ties broken by source: Xu 2004 > Ding 2002 > NBS 514 > PubChem/HSDB.
- Otherwise Gelavizh's `PROPS` value, flagged `props_unverified`: TFP η, Sulfolane ε, 2-Glyme ε,
  3-Glyme, 4-Glyme, DMSO, FEC, MOEMC, DMF ε. These are listed in `smiles/CROSSCHECK_GELAVIZH.md`
  and must be confirmed with her before the paper is written.
- Solvents with neither are excluded from every set: 3-MeSulfolane, DMM, Propylsulfone, Ethyldiglyme,
  Ethylmonoglyme, Freon 11, Pseudocumene. This leaves 31 solvents and 6619 / 6759 train rows, and all
  5959 test rows.

Values are taken at the temperatures their sources report (20–40 °C), not corrected to 25 °C. EC is
at 40 °C (solid at 25 °C) and Sulfolane η at 30 °C (solid at 25 °C).

**Effect on the design.** All 13 primary folds remain. Subset U is unchanged. If a cross-check later
changes a value, the affected results are rerun and reported as post-hoc alongside the originals.
