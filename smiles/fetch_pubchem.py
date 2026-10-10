"""Fetch experimental viscosity and dielectric constant records from PubChem. Run: python3 smiles/fetch_pubchem.py

For each solvent in solvent_properties.csv: SMILES -> PubChem CID, then the PUG-View sections "Viscosity" and
"Dielectric Constant". Every record is kept verbatim with its PubChem reference (which names the original source,
e.g. CRC Handbook or HSDB). Raw JSON goes to smiles/lit_raw/, flat records to smiles/pubchem_records.csv.
No values are chosen here; selection happens by hand in solvent_eta_eps.csv.
"""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
D = HERE.parent / "ca-li-sol-23-challenge"
RAW = HERE / "lit_raw"
API = "https://pubchem.ncbi.nlm.nih.gov/rest"
HEADINGS = ["Viscosity", "Dielectric Constant"]


def get(url):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            time.sleep(2 ** attempt)
        except urllib.error.URLError:
            time.sleep(2 ** attempt)
    raise RuntimeError(f"failed: {url}")


def strings(info):
    v = info.get("Value", {})
    out = [s.get("String", "") for s in v.get("StringWithMarkup", [])]
    if "Number" in v:
        out.append(" ".join(str(x) for x in v["Number"]) + " " + v.get("Unit", ""))
    return [s.strip() for s in out if s.strip()]


def main():
    RAW.mkdir(exist_ok=True)
    sp = pd.read_csv(D / "solvent_properties.csv")
    rows = []
    for r in sp.itertuples():
        q = urllib.parse.quote(r.smiles, safe="")
        js = get(f"{API}/pug/compound/smiles/{q}/cids/JSON")
        cid = js["IdentifierList"]["CID"][0] if js else None
        print(f"{r.name:22s} CID {cid}")
        if not cid:
            continue
        for h in HEADINGS:
            view = get(f"{API}/pug_view/data/compound/{cid}/JSON?heading={urllib.parse.quote(h)}")
            time.sleep(0.25)
            if view is None:
                continue
            (RAW / f"{r.name.replace(' ', '_')}_{h.replace(' ', '_')}.json").write_text(json.dumps(view))
            refs = {x["ReferenceNumber"]: x for x in view["Record"].get("Reference", [])}

            def walk(sec):
                for s in sec.get("Section", []):
                    yield from walk(s)
                if sec.get("TOCHeading") == h:
                    yield from sec.get("Information", [])

            for info in walk(view["Record"]):
                ref = refs.get(info.get("ReferenceNumber"), {})
                for s in strings(info):
                    rows.append({"name": r.name, "cid": cid, "property": h, "value_text": s,
                                 "reference": info.get("Reference", [""])[0] if info.get("Reference") else "",
                                 "source": ref.get("SourceName", ""), "source_url": ref.get("URL", "")})
    out = pd.DataFrame(rows)
    out.to_csv(HERE / "pubchem_records.csv", index=False)
    print(f"\n{len(out)} records; solvents with viscosity: {out[out.property == 'Viscosity'].name.nunique()}, "
          f"with dielectric constant: {out[out.property == 'Dielectric Constant'].name.nunique()} of {len(sp)}")


if __name__ == "__main__":
    main()
