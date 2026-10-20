#!/usr/bin/env bash
# Usage: fetch.sh <url> <sha256> <out>
set -euo pipefail

url=$1
sha256=$2
out=$3

if command -v sha256sum > /dev/null; then
    check() { sha256sum -c -; }
else
    check() { shasum -a 256 -c -; }
fi

mkdir -p "$(dirname "$out")"
trap 'rm -f "$out.tmp"' EXIT
curl -fsSL "$url" -o "$out.tmp"
echo "$sha256  $out.tmp" | check > /dev/null
mv "$out.tmp" "$out"
