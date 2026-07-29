#!/usr/bin/env bash
# Build the UniRef50-specific shared preprocessing cache once before the six-task array.

#$ -l tmem=16G
#$ -l tscratch=120G
#$ -l scratch0free=120G
#$ -l h_rt=72:0:0
#$ -j y
#$ -N homology_u50_cache

set -euo pipefail

FRAMEWORK_ROOT="${FRAMEWORK_SOURCE_ROOT:-${SGE_O_WORKDIR:-$PWD}}"
SAN_ROOT="${SAN_ROOT:-/SAN/bioinf/bmpfp}"
UNIREF50_FASTA="${UNIREF50_FASTA:-$SAN_ROOT/frozen_inputs/uniref50/2026_02/uniref50.fasta.gz}"
IDMAPPING="${IDMAPPING:-$SAN_ROOT/frozen_inputs/idmapping/2026_02/idmapping_selected.tab.gz}"
SPROT="${UNIPROT_SPROT_SEQUENCES:-$SAN_ROOT/frozen_inputs/uniprot/2026_02/uniprot_sprot.dat.gz}"
TREMBL="${UNIPROT_TREMBL_SEQUENCES:-$SAN_ROOT/frozen_inputs/uniprot/2026_02/uniprot_trembl.dat.gz}"
GOA="${GOA:-$SAN_ROOT/frozen_inputs/goa/234/goa_uniprot_all.gaf.234.gz}"
GO_OBO="${GO_OBO:-$SAN_ROOT/frozen_inputs/ontology/2026-06-19/go-basic.obo}"
CACHE_ROOT="${HOMOLOGY_COMMON_PREPROCESSING_CACHE:-$SAN_ROOT/derived_inputs/homology/2026_02/goa_234/sprot-and-trembl/uniref50/common_preprocessing}"
WORK="${WORK_BASE:-/scratch0}/homology_uniref50_common_cache_${JOB_ID:-manual}"

cleanup() {
    local status=$?
    set +e
    cd "$HOME"
    [[ ! -e "$WORK/.homology-uniref50-cache-owned" ]] || rm -rf "$WORK"
    exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT TERM

for path in "$UNIREF50_FASTA" "$IDMAPPING" "$SPROT" "$TREMBL" "$GOA" "$GO_OBO"; do
    [[ -s "$path" ]] || { echo "Required frozen input is missing or empty: $path" >&2; exit 1; }
done
[[ -d "$SAN_ROOT" && -w "$SAN_ROOT" ]] || {
    echo "SAN root is unavailable or not writable: $SAN_ROOT" >&2
    exit 1
}
[[ ! -e "$WORK" ]] || { echo "Refusing to reuse scratch directory: $WORK" >&2; exit 1; }
mkdir -p "$WORK/contracts" "$WORK/build"
touch "$WORK/.homology-uniref50-cache-owned"

# shellcheck source=../../scripts/reproduction_common.sh
source "$FRAMEWORK_ROOT/scripts/reproduction_common.sh"
load_framework_paths "$FRAMEWORK_ROOT"
add_mmfp_singularity_bind "$SAN_ROOT"
activate_or_create_mmfp_env
PYTHON_BIN="$(command -v python)"
export PYTHONPATH="$FRAMEWORK_ROOT/benchmark_builders/homology_cluster/src${PYTHONPATH:+:$PYTHONPATH}"
FRAMEWORK_REVISION="$(cd "$FRAMEWORK_ROOT" && git rev-parse HEAD)"

if [[ -s "$CACHE_ROOT/CACHE_COMPLETE.json" ]] && \
   "$PYTHON_BIN" -m homology_cluster_benchmark.common_cache verify \
       --cache-dir "$CACHE_ROOT" \
       --source-scope sprot-and-trembl \
       --uniref-level 50 \
       --full-hashes; then
    echo "Validated UniRef50 common preprocessing cache already exists: $CACHE_ROOT"
    exit 0
fi

MANIFEST="$WORK/contracts/frozen_input_manifest.json"
POLICY="$WORK/contracts/unused_runtime_policy.json"
"$PYTHON_BIN" -m homology_cluster_benchmark.runtime_contract prepare \
    --manifest-out "$MANIFEST" \
    --policy-out "$POLICY" \
    --source-scope sprot-and-trembl \
    --framework-revision "$FRAMEWORK_REVISION" \
    --uniref-level 50 \
    --uniref50-fasta "$UNIREF50_FASTA" \
    --uniref50-fasta-url https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta.gz \
    --uniref50-fasta-acquisition provided-persistent-store \
    --idmapping "$IDMAPPING" \
    --idmapping-url https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz \
    --idmapping-acquisition provided-persistent-store \
    --uniprot-sprot-sequences "$SPROT" \
    --uniprot-sprot-sequences-url https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz \
    --uniprot-sprot-sequences-acquisition provided-persistent-store \
    --uniprot-trembl-sequences "$TREMBL" \
    --uniprot-trembl-sequences-url https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz \
    --uniprot-trembl-sequences-acquisition provided-persistent-store \
    --goa "$GOA" \
    --goa-url https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/goa_uniprot_all.gaf.234.gz \
    --goa-acquisition provided-persistent-store \
    --go-obo "$GO_OBO" \
    --go-obo-url https://release.geneontology.org/2026-06-19/ontology/go-basic.obo \
    --go-obo-acquisition provided-persistent-store

"$PYTHON_BIN" -m homology_cluster_benchmark.common_cache build \
    --output-dir "$CACHE_ROOT" \
    --work-dir "$WORK/build" \
    --frozen-input-manifest "$MANIFEST" \
    --source-scope sprot-and-trembl \
    --uniref-level 50 \
    --uniref50-fasta "$UNIREF50_FASTA" \
    --idmapping "$IDMAPPING" \
    --uniprot-sprot-sequences "$SPROT" \
    --uniprot-trembl-sequences "$TREMBL" \
    --goa "$GOA" \
    --go-obo "$GO_OBO" \
    --replace-existing

"$PYTHON_BIN" -m homology_cluster_benchmark.common_cache verify \
    --cache-dir "$CACHE_ROOT" \
    --source-scope sprot-and-trembl \
    --uniref-level 50 \
    --full-hashes
echo "Published validated UniRef50 common preprocessing cache: $CACHE_ROOT"
