"""Build the frozen set-A property table. Run: python3 smiles/build_props.py

M and rho: from the competition files (smiles/solvent_props_kaggle.csv, Amendment 1).
eps and eta: a cited literature value from smiles/solvent_eta_eps_lit.csv where one exists, chosen as the value
measured closest to 25 C, ties broken by source (Xu 2004 > Ding 2002 > NBS 514 > PubChem/HSDB); otherwise
Gelavizh's PROPS value (pipeline.py), flagged "props_unverified" for cross-checking with her (Amendment 2). Solvents with neither are left empty and their rows
are excluded from every feature set.
Writes smiles/solvent_props_final.csv.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent))
from pipeline import PROPS  # noqa: E402

PRIORITY = ["Xu 2004", "Ding Xu Jow 2002", "Maryott & Smith 1951", "PubChem"]
PROPS_IDX = {"eps": 2, "eta": 3}


def rank(source):
    return next((i for i, p in enumerate(PRIORITY) if source.startswith(p)), len(PRIORITY))


def main():
    kag = pd.read_csv(HERE / "solvent_props_kaggle.csv")
    lit = pd.read_csv(HERE / "solvent_eta_eps_lit.csv")
    lit["rank"] = lit.source.map(rank)
    lit["dT"] = (lit.T_C - 25).abs()
    rows = []
    for name in kag.name:
        r = {"name": name, "M_g_mol": kag.set_index("name").M_g_mol[name], "rho_g_cm3": kag.set_index("name").rho_g_cm3[name]}
        for prop in ["eps", "eta"]:
            cand = lit[(lit.name == name) & (lit.property == prop)].sort_values(["dT", "rank"])
            props_val = PROPS[name][PROPS_IDX[prop]] if name in PROPS else np.nan
            if len(cand):
                c = cand.iloc[0]
                r.update({prop: c.value, f"{prop}_T_C": c.T_C, f"{prop}_source": f"{c.source}; {c.locator}",
                          f"{prop}_status": "verified"})
            elif not np.isnan(props_val):
                r.update({prop: props_val, f"{prop}_T_C": 25 if name != "EC" else 40,
                          f"{prop}_source": "pipeline.py PROPS (Gelavizh)", f"{prop}_status": "props_unverified"})
            else:
                r.update({prop: np.nan, f"{prop}_T_C": np.nan, f"{prop}_source": "", f"{prop}_status": "missing"})
            r[f"{prop}_props"] = props_val
        rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(HERE / "solvent_props_final.csv", index=False)

    tr = pd.read_csv(HERE.parent / "ca-li-sol-23-challenge" / "train.csv")
    usable = set(out.name[out.eps.notna() & out.eta.notna()])
    covered = tr.solvents.str.split(";").apply(lambda l: all(s in usable for s in l))
    pd.set_option("display.width", 200)
    print(out[["name", "eps", "eps_status", "eta", "eta_status"]].to_string(index=False))
    print(f"\nusable solvents {len(usable)} / {len(out)}; train rows covered {covered.sum()} / {len(tr)}")
    for prop in ["eps", "eta"]:
        d = out[(out[f"{prop}_status"] == "verified") & out[f"{prop}_props"].notna()]
        rel = (d[prop] / d[f"{prop}_props"] - 1).abs()
        print(f"{prop}: verified vs PROPS differ by > 3%:", dict(zip(d.name[rel > 0.03], d[prop][rel > 0.03])))


if __name__ == "__main__":
    main()
