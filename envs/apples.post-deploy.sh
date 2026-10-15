#!/usr/bin/env bash
# Run once after `micromamba create` from envs/apples.linux-64.pin.txt, with the env activated. Idempotent.
set -euo pipefail

pip install --no-deps --require-hashes --no-build-isolation -r /dev/stdin <<'EOF'
apples==2.0.11 --hash=sha256:0df2367ceac167c6c347954911643a0305c2f293a5f31a76bc86bd568bf12115
EOF
