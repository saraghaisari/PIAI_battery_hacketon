"""Pure-solvent molar mass and density, recovered from the competition files alone. Run: python3 smiles/props_kaggle.py

M: RDKit average molecular weight from solvent_properties.csv SMILES.
rho: metaData.csv says mixture_density_g_cm3 is the ideal additive-volume density of the pure solvents,
     1/rho_mix = sum_i w_i / rho_i (w = mass fractions). Every distinct mixture in train+test gives one linear
     equation in the 38 unknowns 1/rho_i; solved by least squares. Uses no target values.
Writes smiles/solvent_props_kaggle.csv and compares with Gelavizh's PROPS (pipeline.py).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / "ca-li-sol-23-challenge"
sys.path.insert(0, str(ROOT))
from pipeline import PROPS  # noqa: E402


def main():
    sp = pd.read_csv(D / "solvent_properties.csv")
    df = pd.concat([pd.read_csv(D / "train.csv"), pd.read_csv(D / "test.csv")])
    M = {r.name: Descriptors.MolWt(Chem.MolFromSmiles(r.smiles)) for r in sp.itertuples()}
    names = sorted(M)
    idx = {n: i for i, n in enumerate(names)}

    mix = df[["solvents", "solvent_fracs_mol", "mixture_density_g_cm3"]].drop_duplicates()
    A, b = np.zeros((len(mix), len(names))), np.zeros(len(mix))
    for j, r in enumerate(mix.itertuples()):
        sv = r.solvents.split(";")
        w = np.array(r.solvent_fracs_mol.split(";"), float) * np.array([M[n] for n in sv])
        for n, wi in zip(sv, w / w.sum()):
            A[j, idx[n]] += wi
        b[j] = 1.0 / r.mixture_density_g_cm3
    inv_rho, _, rank, _ = np.linalg.lstsq(A, b, rcond=None)
    resid = np.abs(A @ inv_rho - b).max()
    print(f"{len(mix)} distinct mixtures, rank {rank} of {len(names)}, max |residual| {resid:.1e}")
    assert rank == len(names) and resid < 1e-6, "densities not uniquely recovered"

    out = pd.DataFrame({"name": names, "M_g_mol": [M[n] for n in names], "rho_g_cm3": 1.0 / inv_rho})
    out["M_props"] = [PROPS[n][0] if n in PROPS else np.nan for n in names]
    out["rho_props"] = [PROPS[n][1] if n in PROPS else np.nan for n in names]
    out.round(4).to_csv(Path(__file__).parent / "solvent_props_kaggle.csv", index=False)
    pd.set_option("display.width", 200)
    print(out.round(3).to_string(index=False))
    d = (out.rho_g_cm3 - out.rho_props).abs()
    print("\nPROPS density differs by > 0.02 g/cm3:", out.name[d > 0.02].tolist())


if __name__ == "__main__":
    main()
