#!/bin/bash
# verify_cafa3_csvs.sh — MAXIMALLY thorough provenance check:
# does Zijan's mmfp_data_splits regenerate EXACTLY from Zenodo record 7409660?
# Fully isolated in VERIFY_CSV_WORKDIR. PFP is only READ, never written.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
if [ -f "${HERE}/configs/paths.local.sh" ]; then
  # Machine-specific paths are intentionally not committed.
  # shellcheck disable=SC1091
  source "${HERE}/configs/paths.local.sh"
fi

VDIR="${VERIFY_CSV_WORKDIR:-${HOME}/mmfp_csv_verify}"
RAW="${VDIR}/cafa3_raw"
GEN="${VDIR}/generated_splits"
PFP="${PFP_DIR:-${HOME}/PFP}"
MMFP_ENV="${MMFP_ENV:-mmfp}"

# --- Activate the env FIRST, before any python call ---
eval "$(micromamba shell hook --shell bash)"
micromamba activate "${MMFP_ENV}"
python --version   # sanity: should print the mmfp env's Python, not pyenv

mkdir -p "${RAW}" "${GEN}"
echo "==> Verification sandbox: ${VDIR}  (PFP is read-only here)"
echo "==> PFP reference path: ${PFP}"
echo "==> MMFP environment: ${MMFP_ENV}"

# --- 1. Download the 9 required CSVs (skip augmented; not used by the script) ---
cd "${RAW}"
BASE="https://zenodo.org/records/7409660/files"
cat > "${VDIR}/zenodo_md5.txt" <<'EOF'
e9a4b239cd47a7ac80975f63e259581e  bp-test.csv
85c19594547a503956226b9c225efc5d  bp-training.csv
c2674223770d6a8cf680dd9335d51ebe  bp-validation.csv
0e5dc8528ca95e8897b10cddaa12a775  cc-test.csv
074b13dd50fad4a6a4f13e4d8d4105d6  cc-training.csv
cdc8ceefcab4fb8c9278dd07c184327f  cc-validation.csv
2735e408dd57f6de29b1538f6b150d68  mf-test.csv
b31a8f22b5934aef61b76ec3b89296da  mf-training.csv
897921ce5df8174672200320926ccc87  mf-validation.csv
EOF
for aspect in bp cc mf; do
  for split in training validation test; do
    f="${aspect}-${split}.csv"
    [ -f "$f" ] || { echo "==> Downloading $f"; wget -c "${BASE}/${f}?download=1" -O "$f"; }
  done
done
echo "==> Verifying downloaded CSV md5s against Zenodo..."
if command -v md5sum >/dev/null 2>&1; then
  md5sum -c "${VDIR}/zenodo_md5.txt" || { echo "CSV md5 MISMATCH — download corrupt or record changed."; exit 1; }
else
  while read -r want name; do
    got=$(md5 -q "$name"); [ "$got" = "$want" ] && echo "  OK  $name" || { echo "  BAD $name ($got != $want)"; exit 1; }
  done < "${VDIR}/zenodo_md5.txt"
fi
echo "==> All 9 CSVs authenticated."

# --- 2. Copy the prep script into the sandbox AND patch the column-name quirk ---
#     The MF CSVs use 'protein' (singular); BP/CC use 'proteins'. The published
#     script hard-codes 'proteins' and crashes on MF. Normalise after read_csv.
if [ ! -f "${PFP}/scripts/prepare_cafa3_data.py" ]; then
  echo "Missing PFP prepare script: ${PFP}/scripts/prepare_cafa3_data.py" >&2
  echo "Set PFP_DIR in configs/paths.local.sh or the environment." >&2
  exit 1
fi
cp "${PFP}/scripts/prepare_cafa3_data.py" "${VDIR}/prepare_cafa3_data.py"
python - "${VDIR}/prepare_cafa3_data.py" <<'PY'
import sys
from pathlib import Path
p = Path(sys.argv[1]); src = p.read_text()
needle = "        df = pd.read_csv(filepath)\n"
inject = needle + "        df = df.rename(columns={'protein': 'proteins'})  # normalise MF singular -> plural\n"
assert needle in src, "read_csv line not found — script changed upstream?"
p.write_text(src.replace(needle, inject, 1))
print("==> Patched sandbox prepare_cafa3_data.py: normalised 'protein' -> 'proteins'.")
PY

cd "${VDIR}"
python prepare_cafa3_data.py --cafa3-dir "${RAW}" --output-dir "${GEN}"

# --- 3. EXHAUSTIVE comparison against Zijan's reference (PFP/data, read-only) ---
python - "$GEN" "${PFP}/data" <<'PY'
import sys, json, numpy as np, scipy.sparse as ssp
from pathlib import Path
gen, ref = Path(sys.argv[1]), Path(sys.argv[2])
aspects, splits = ["BPO","CCO","MFO"], ["train","valid","test"]
results = []
def rec(label, ok, detail=""): results.append((label, ok, detail))
def status(ok): return "MATCH" if ok else "DIFFER"
for a in aspects:
    gt_g, gt_r = gen/f"{a}_go_terms.json", ref/f"{a}_go_terms.json"
    if gt_g.exists() and gt_r.exists():
        Lg, Lr = json.load(open(gt_g)), json.load(open(gt_r))
        rec(f"{a} go_terms (exact order)", Lg == Lr,
            f"gen={len(Lg)} ref={len(Lr)} set_equal={set(Lg)==set(Lr)}")
    else:
        rec(f"{a} go_terms", False, "file missing")
    for s in splits:
        ng, nr = gen/f"{a}_{s}_names.npy", ref/f"{a}_{s}_names.npy"
        if ng.exists() and nr.exists():
            g = np.load(ng, allow_pickle=True); r = np.load(nr, allow_pickle=True)
            exact = np.array_equal(g, r); set_eq = set(g.tolist()) == set(r.tolist())
            rec(f"{a}_{s} names (exact+order)", exact,
                f"gen={len(g)} ref={len(r)} set_equal={set_eq} order_equal={exact}")
        else:
            rec(f"{a}_{s} names", False, "file missing"); continue
        sg, sr = gen/f"{a}_{s}_sequences.json", ref/f"{a}_{s}_sequences.json"
        if sg.exists() and sr.exists():
            Dg, Dr = json.load(open(sg)), json.load(open(sr))
            same_keys = set(Dg)==set(Dr); same_vals = same_keys and all(Dg[k]==Dr[k] for k in Dg)
            rec(f"{a}_{s} sequences", same_vals, f"keys_equal={same_keys} values_equal={same_vals}")
        else:
            rec(f"{a}_{s} sequences", False, "file missing")
        lg, lr = gen/f"{a}_{s}_labels.npz", ref/f"{a}_{s}_labels.npz"
        if lg.exists() and lr.exists():
            Mg, Mr = ssp.load_npz(lg), ssp.load_npz(lr)
            if Mg.shape != Mr.shape:
                rec(f"{a}_{s} labels", False, f"shape gen={Mg.shape} ref={Mr.shape}")
            else:
                diff = (Mg != Mr).nnz
                rec(f"{a}_{s} labels (exact)", diff==0, f"shape={Mg.shape} differing_entries={diff}")
        else:
            rec(f"{a}_{s} labels", False, "file missing")
print("\n================ EXHAUSTIVE VERIFICATION ================")
width = max(len(l) for l,_,_ in results)
for label, ok, detail in results:
    print(f"  [{status(ok):6}] {label:<{width}}  {detail}")
allok = all(ok for _,ok,_ in results)
print("\n==> OVERALL:",
      "FULL BIT-LEVEL MATCH — Zenodo 7409660 is conclusively Zijan's source."
      if allok else
      "NOT a full match — see DIFFER rows above for exactly which artifact and how.")
PY

echo ""
echo "==> Done. All artifacts in ${VDIR}. PFP was never modified."
