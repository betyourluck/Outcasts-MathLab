#!/bin/sh
# Regenerate every output of the failure atlas from scratch (Python 3 standard library only).
set -e
cd "$(dirname "$0")"
python failure_atlas.py > summary.json
gzip -n -9 -f failure_atlas_rows.json
python h2prime.py > h2prime_summary.json
gzip -n -9 -f h2prime_rows.json
python exponent_atlas.py > exponent_summary.json
gzip -n -9 -f exponent_rows.json
python h10prime.py > h10prime_summary.json
python h11_independence.py > h11_summary.json
python h11_posthoc.py > h11_posthoc.json
# auxiliary checks (printed only)
gzip -dc failure_atlas_rows.json.gz > failure_atlas_rows.json
python h4_chars.py
rm failure_atlas_rows.json
python jacobi_check.py
