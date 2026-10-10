# To cross-check with Gelavizh (IMPORTANT, open)

Set A (physics) currently uses her `pipeline.py` `PROPS` values wherever we found no cited source
(Amendment 2). Each item below must be confirmed or replaced, with a source, before the paper is written.
Any change after training starts is reported as a post-hoc change.

## Her values used without a source (status `props_unverified` in `solvent_props_final.csv`)

| Solvent | Value used | Train rows | Own validation fold | Question |
|---|---|---|---|---|
| TFP | η 2.5 cP | 540 | yes | Where does 2.5 cP come from? Ding 2002 and Xu 2004 do not report it. |
| Sulfolane | ε 43.3 | 334 | yes | Source? (η 10.34 cP at 30 °C is verified) |
| 3-Glyme | ε 7.5, η 1.96 | 217 | yes | Source? |
| 2-Glyme | ε 7.4 | 145 | yes | Source? (η verified, 1.089 cP at 20 °C) |
| DMSO | ε 46.5, η 1.99 | 76 | yes | Source? PubChem has 2.47 cP at 20 °C (HSDB) vs 1.996 cP (Wikipedia) |
| FEC | ε 78.4, η 4.1 | 42 | no | Source? Her comment marks it least certain |
| MOEMC | ε 9.0, η 1.5 | 24 | no | Source? Her comment marks it least certain |
| 4-Glyme | ε 7.8, η 3.39 | 36 | no | Source? |
| DMF | ε 36.7 | 6 | no | Source? (η verified) |

## Her values that disagree with a cited source (we use the cited value)

| Solvent | Hers | Cited | Source |
|---|---|---|---|
| TFP | ε 11.0 | ε 10.25 at 25 °C (10.5 at 20 °C) | Ding, Xu & Jow 2002; Xu 2004 Table 14 |
| TFP | ρ 1.59 | ρ 1.487 | recovered from Kaggle mixture densities |
| MOEMC | ρ 1.10 | ρ 1.070 | recovered from Kaggle mixture densities |
| 2-MeTHF | ε 6.97 | ε 6.2 | Xu 2004 Table 2 |
| EA | η 0.43 | η 0.45 | Xu 2004 Table 1 (HSDB: 0.423) |
| 2-Glyme | η 0.99 at 25 °C | η 1.089 at 20 °C | PubChem/HSDB (different temperature) |
| AN | ε 35.9 | ε 37.5 at 20 °C | NBS 514 (different temperature) |

## Also ask

- `pip freeze` from her runs (scikit-learn version explains different baseline numbers).
- Kaggle leaderboard score of `submission_final.csv`.
- Access to CRC Handbook / Riddick for the open items above.
