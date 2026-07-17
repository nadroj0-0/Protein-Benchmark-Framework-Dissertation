#!/usr/bin/env bash
# Portable lookup helpers for large, pre-downloaded research artifacts.
#
# Resolution precedence is deliberately simple:
#   1. an explicit caller-supplied path, when it exists;
#   2. the matching path in ARTIFACT_CATALOG, when it exists;
#   3. no result, so the calling workflow can use its original download path.
#
# Catalogue files are parsed as data and are never sourced as shell code.

artifact_catalog_warn() {
  printf 'Artifact catalogue warning: %s\n' "$*" >&2
}

artifact_catalog_configure() {
  local framework_root="$1"
  local requested="${2:-${ARTIFACT_CATALOG:-}}"
  local local_catalog="${framework_root}/configs/artifact_paths.local.tsv"

  if [[ -z "$requested" && -f "$local_catalog" ]]; then
    requested="$local_catalog"
  fi
  if [[ -z "$requested" ]]; then
    ARTIFACT_CATALOG=""
    export ARTIFACT_CATALOG
    return 0
  fi
  if [[ ! -f "$requested" ]]; then
    artifact_catalog_warn "catalogue does not exist; download fallbacks remain enabled: $requested"
    ARTIFACT_CATALOG=""
    export ARTIFACT_CATALOG
    return 0
  fi

  ARTIFACT_CATALOG="$(cd "$(dirname "$requested")" && pwd -P)/$(basename "$requested")"
  artifact_catalog_validate "$ARTIFACT_CATALOG" || return 1
  export ARTIFACT_CATALOG
}

artifact_catalog_validate() {
  local catalog="$1"
  awk -F '\t' '
    BEGIN { valid=1 }
    /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
    $1 == "artifact_id" && $2 == "path" { next }
    NF < 2 || $1 == "" || $2 == "" {
      printf "Malformed artifact catalogue row %d in %s\n", NR, FILENAME > "/dev/stderr"
      valid=0
      next
    }
    $1 !~ /^[A-Za-z0-9._-]+$/ {
      printf "Unsafe artifact ID on row %d in %s: %s\n", NR, FILENAME, $1 > "/dev/stderr"
      valid=0
    }
    $2 !~ /^\// {
      printf "Artifact path is not absolute on row %d in %s: %s\n", NR, FILENAME, $2 > "/dev/stderr"
      valid=0
    }
    seen[$1]++ {
      printf "Duplicate artifact ID on row %d in %s: %s\n", NR, FILENAME, $1 > "/dev/stderr"
      valid=0
    }
    END { exit !valid }
  ' "$catalog"
}

artifact_catalog_lookup() {
  local artifact_id="$1"
  [[ -n "${ARTIFACT_CATALOG:-}" && -f "$ARTIFACT_CATALOG" ]] || return 1
  awk -F '\t' -v wanted="$artifact_id" '
    /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
    $1 == wanted { print $2; found=1; exit }
    END { exit !found }
  ' "$ARTIFACT_CATALOG"
}

resolve_artifact_path() {
  local artifact_id="$1"
  local explicit_path="${2:-}"
  local catalog_path=""

  if [[ -n "$explicit_path" ]]; then
    if [[ -s "$explicit_path" ]]; then
      printf '%s\n' "$explicit_path"
      return 0
    fi
    artifact_catalog_warn "explicit path for $artifact_id is missing or empty; trying catalogue/download fallback: $explicit_path"
  fi

  catalog_path="$(artifact_catalog_lookup "$artifact_id" 2>/dev/null || true)"
  if [[ -n "$catalog_path" ]]; then
    if [[ -s "$catalog_path" ]]; then
      printf '%s\n' "$catalog_path"
      return 0
    fi
    artifact_catalog_warn "catalogue path for $artifact_id is missing or empty; download fallback remains enabled: $catalog_path"
  fi
  return 1
}

artifact_catalog_bind_parent() {
  local artifact_id="$1"
  local explicit_path="${2:-}"
  local resolved=""
  resolved="$(resolve_artifact_path "$artifact_id" "$explicit_path" 2>/dev/null || true)"
  [[ -n "$resolved" ]] || return 0
  add_mmfp_singularity_bind "$(dirname "$resolved")"
}

canonical_cafa3_artifact_id() {
  local name="$1"
  name="${name%.csv}"
  name="${name//-/_}"
  case "$name" in
    bp_training|bp_validation|bp_test|cc_training|cc_validation|cc_test|\
    mf_training|mf_validation|mf_test)
      printf 'cafa3_%s\n' "$name"
      ;;
    *) return 1 ;;
  esac
}

zijian_embedding_artifact_id() {
  case "$1" in
    mmfp_embeddings_prott5.tar.gz) printf '%s\n' zijian_prott5_embeddings ;;
    mmfp_embeddings_struct_ppi.tar.gz) printf '%s\n' zijian_structure_ppi_embeddings ;;
    mmfp_embeddings_text_temporal.tar.gz) printf '%s\n' zijian_text_embeddings ;;
    *) return 1 ;;
  esac
}

zijian_bundle_artifact_id() {
  case "$1" in
    mmfp_embeddings_prott5|mmfp_embeddings_prott5.tar.gz)
      printf '%s\n' zijian_prott5_embeddings ;;
    mmfp_embeddings_struct_ppi|mmfp_embeddings_struct_ppi.tar.gz)
      printf '%s\n' zijian_structure_ppi_embeddings ;;
    mmfp_embeddings_text_temporal|mmfp_embeddings_text_temporal.tar.gz)
      printf '%s\n' zijian_text_embeddings ;;
    mmfp_checkpoints|mmfp_checkpoints.tar.gz)
      printf '%s\n' zijian_checkpoints ;;
    mmfp_data_splits|mmfp_data_splits.tar.gz)
      printf '%s\n' zijian_data_splits ;;
    *) return 1 ;;
  esac
}

stage_or_download_artifact() {
  local artifact_id="$1"
  local explicit_path="$2"
  local destination="$3"
  local url="$4"
  local source_path=""
  source_path="$(resolve_artifact_path "$artifact_id" "$explicit_path" || true)"
  if [[ -n "$source_path" ]]; then
    echo "Staging existing artifact: $source_path"
    cp -p "$source_path" "$destination"
    return 0
  fi
  echo "Downloading artifact: $url"
  if command -v wget >/dev/null 2>&1; then
    wget -c "$url" -O "$destination"
  elif command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 5 --continue-at - --output "$destination" "$url"
  else
    echo "Neither wget nor curl is available" >&2
    return 1
  fi
}
