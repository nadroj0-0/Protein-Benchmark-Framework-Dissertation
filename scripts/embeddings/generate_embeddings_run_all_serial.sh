#!/bin/bash
# generate_embeddings_run_all.sh  — SUB-ORCHESTRATOR (not a standalone entry point)
# Called by reproduce_embeddings_retrain_eval.sh AFTER it has cloned PFP, built the
# env, and cd'd into the repo root. Assumes CWD is the PFP repo root and the mmfp
# env is active. Downloads external deps, regenerates splits, then generates all 4
# modality embeddings in serial. (Modalities are separate scripts so they can be
# parallelised later.)
#
# Builds data/proteins.fasta before the ProtT5 stage so FASTA IDs match split IDs.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../.." && pwd)"
DEPENDENCY_ENV="${DEPENDENCY_ENV:-external/dependency_env.sh}"

# --- External database dependencies (writes external/dependency_env.sh) ---
echo "==> [0/8] External dependencies"
bash "${HERE}/generate_embeddings_dependencies.sh"
if [ ! -f "${DEPENDENCY_ENV}" ]; then
    echo "Missing dependency environment file: ${DEPENDENCY_ENV}" >&2
    echo "Expected generate_embeddings_dependencies.sh to create it from the PFP repo root." >&2
    exit 1
fi
# shellcheck disable=SC1090
source "${DEPENDENCY_ENV}"
echo "==> Loaded dependency environment: ${DEPENDENCY_ENV}"

echo "==> [1/8] Data preparation (raw CSVs -> splits/labels/sequences)"
bash "${HERE}/generate_embeddings_prepare_data.sh"

echo "[2/8] Verifying generated dataset splits..."
python "${REPO_ROOT}/scripts/verification/verify_splits.py" \
    --data-dir data \
    --strict

echo "==> [3/8] Building proteins.fasta from split sequences"
python "${REPO_ROOT}/scripts/embeddings/generate_embeddings_fasta.py" --data-dir data --config "${REPO_ROOT}/configs/cafa3.json"

echo "==> [4/8] Sequence (ProtT5) embeddings"
bash "${HERE}/generate_embeddings_sequence.sh"

echo "==> [5/8] Text embeddings"
bash "${HERE}/generate_embeddings_text.sh"

echo "==> [6/8] Structure embeddings"
bash "${HERE}/generate_embeddings_structure.sh"

echo "==> [7/8] PPI embeddings"
bash "${HERE}/generate_embeddings_ppi.sh"

## --- [4-7/8] Parallel modality embeddings ---------------------------------
## prott5 / text / structure -> one GPU each;  ppi -> CPU, concurrent.
#
#LOGDIR="results/embedding_logs"
#mkdir -p "$LOGDIR"
#
## SGE's gpu PE exposes the allocated devices via CUDA_VISIBLE_DEVICES.
## Split them; fall back to 0,1,2 for interactive testing.
#IFS=',' read -ra GPUS <<< "${CUDA_VISIBLE_DEVICES:-0,1,2}"
#gpu() { echo "${GPUS[$1]:-${GPUS[0]}}"; }
#echo "==> Parallel embeddings on GPUs: ${GPUS[*]} (ppi on CPU)"
#
#CUDA_VISIBLE_DEVICES="$(gpu 0)" DEVICE=cuda \
#    bash "${HERE}/generate_embeddings_sequence.sh"  > "$LOGDIR/seq.log"    2>&1 &
#PID_SEQ=$!
#CUDA_VISIBLE_DEVICES="$(gpu 1)" \
#    bash "${HERE}/generate_embeddings_text.sh"      > "$LOGDIR/text.log"   2>&1 &
#PID_TXT=$!
#CUDA_VISIBLE_DEVICES="$(gpu 2)" DEVICE=cuda \
#    bash "${HERE}/generate_embeddings_structure.sh" > "$LOGDIR/struct.log" 2>&1 &
#PID_STR=$!
#CUDA_VISIBLE_DEVICES="" \
#    bash "${HERE}/generate_embeddings_ppi.sh"       > "$LOGDIR/ppi.log"    2>&1 &
#PID_PPI=$!
#
## Join, collecting each exit status. Using `wait` as an `if` condition
## suppresses `set -e`, so one failure doesn't abort the other waits.
#rc=0
#for np in "sequence:$PID_SEQ" "text:$PID_TXT" "structure:$PID_STR" "ppi:$PID_PPI"; do
#    name="${np%%:*}"; pid="${np##*:}"
#    if wait "$pid"; then echo "==> [$name] OK"
#    else echo "==> [$name] FAILED (logs/${name}.log)"; rc=1; fi
#done
#[ "$rc" -eq 0 ] || { echo "A modality job failed; aborting before verify."; exit 1; }


echo "[8/8] Verifying generated embeddings..."
python "${REPO_ROOT}/scripts/verification/verify_embeddings.py" \
    --data-dir data \
    --config "${REPO_ROOT}/configs/cafa3.json" \
    --strict

echo "==> Embedding generation complete. data/embedding_cache/ populated."
