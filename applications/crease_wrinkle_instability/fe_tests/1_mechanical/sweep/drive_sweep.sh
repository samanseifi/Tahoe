#!/usr/bin/env bash
cd "$(dirname "$0")"
{
echo "0.0  0.56"
echo "0.25 0.65"
echo "0.5  0.69"
echo "1.0  0.75"
echo "2.0  0.81"
} | xargs -P 5 -L1 bash -c 'python3 sweep.py "$@"' _
echo "=== SWEEP DONE ==="
python3 extract_sweep.py > sweep_results.txt 2>&1
echo "=== FINALIZED ===" >> sweep_results.txt
