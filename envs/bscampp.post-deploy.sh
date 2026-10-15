#!/usr/bin/env bash
# Run once after `micromamba create` from envs/bscampp.linux-64.pin.txt, with the env activated. Idempotent.
set -euo pipefail

pip install --no-deps --require-hashes -r /dev/stdin <<'EOF'
bscampp==1.0.8 --hash=sha256:825aebbe464bae36cef448f6a94740119bb04946167260eae847cbe7d58c2389
EOF

SP=$(python -c 'import bscampp,os;print(os.path.dirname(bscampp.__file__))')
CONFIG_DIR="$CONDA_PREFIX/share/bscampp"

mkdir -p "$CONFIG_DIR"
cat > "$CONFIG_DIR/main.config" <<EOF
[basic]
pplacer_path = $CONDA_PREFIX/bin/pplacer
epang_path = $CONDA_PREFIX/bin/epa-ng
taxit_path = $CONDA_PREFIX/bin/taxit
hamming_distance_dir = $SP/tools/hamming_distance
EOF
echo "$CONFIG_DIR" > "$SP/home.path"

# bscampp functions.py ensureBinaryExecutable treats exit 255 as "works"; anything else
# makes it recompile the bundled binaries into site-packages at runtime.
for b in fragment_hamming hamming homology; do
    rc=0
    "$SP/tools/hamming_distance/$b" >/dev/null 2>&1 || rc=$?
    if [ "$rc" -ne 255 ]; then
        echo "bundled $b exited $rc, expected 255; bscampp would recompile at runtime" >&2
        exit 1
    fi
done
