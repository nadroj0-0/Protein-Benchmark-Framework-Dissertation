#!/usr/bin/env bash

VERSIONS="0.41.2 0.41.1 0.40.0 0.39.0 0.38.0 0.37.0 0.36.1 0.36.0 0.35.0 0.34.1 0.34.0 0.33.0 0.32.0 0.31.0 0.30.0 0.29.0 0.28.0 0.27.0 0.26.0 0.25.0 0.24.0 0.23.0 0.22.0 0.21.0 0.20.1 0.20.0"

echo "Testing Biotite versions..."

for VER in $VERSIONS
do
    echo "===================================================="
    echo "Testing biotite==$VER"
    echo "===================================================="

    pip uninstall -y biotite biotraj >/dev/null 2>&1

    if ! pip install --only-binary=:all: "biotite==$VER" >/tmp/biotite_install.log 2>&1
    then
        echo "INSTALL FAILED"
        continue
    fi

    python - <<'PY'
import sys

try:
    from biotite.structure import filter_backbone
except Exception as e:
    print("FAIL filter_backbone:", e)
    sys.exit(1)

try:
    import esm.inverse_folding
except Exception as e:
    print("FAIL esm.inverse_folding:", e)
    sys.exit(2)

try:
    import esm
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
except Exception as e:
    print("FAIL model:", e)
    sys.exit(3)

print("SUCCESS")
PY

    if [ "$?" -eq 0 ]; then
        echo
        echo "FOUND WORKING VERSION: biotite==$VER"
        exit 0
    fi
done

echo "No version worked."
exit 1
