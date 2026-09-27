#!/usr/bin/env bash
# run_smoke.sh — run one tahoe command (serial or under mpirun) and fail if it crashes OR aborts.
# tahoe exits 0 even when a run gives up (issue #74), so the exit code alone is not enough: the
# output is also checked for "Exiting time sequence" (load-step cutting exhausted) and
# "exit on exception". Usage: run_smoke.sh <command...>, e.g.
#   run_smoke.sh ./build/bin/tahoe -f deck.xml
#   run_smoke.sh mpirun -np 2 ./build/bin/tahoe -f deck.xml
log=$(mktemp)
"$@" 2>&1 | tee "$log"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]; then
    echo "run_smoke.sh: command exited with code $rc"; rm -f "$log"; exit "$rc"
fi
if grep -q "Exiting time sequence\|exit on exception" "$log"; then
    echo "run_smoke.sh: the run aborted (see 'Exiting time sequence' / 'exit on exception' above)"
    rm -f "$log"; exit 1
fi
rm -f "$log"
