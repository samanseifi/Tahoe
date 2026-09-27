#!/usr/bin/env bash
# run_benchmarks.sh — Run all Tahoe benchmark levels and report failures
# Usage: ./run_benchmarks.sh [level.0] [level.1] ...  (default: all levels)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TAHOE="${TAHOE:-$SCRIPT_DIR/build/bin/tahoe}"
COMPARE="${COMPARE:-$SCRIPT_DIR/build/bin/compare}"
BENCH_ROOT="${BENCH_ROOT:-$SCRIPT_DIR/benchmark_XML}"
MPI_NP="${MPI_NP:-4}"        # number of MPI ranks for level.3; override via env
MPIFLAGS="${MPIFLAGS:-}"     # extra mpirun flags, e.g. "--oversubscribe" on CI

# Detect system MPI (avoid conda wrappers which use a different MPI implementation)
# Override via MPIRUN env variable.
if [ -z "${MPIRUN+x}" ]; then
    MPIRUN=""
    for candidate in /usr/bin/mpirun /usr/local/bin/mpirun /usr/bin/mpirun.openmpi; do
        if [ -x "$candidate" ]; then
            MPIRUN="$candidate"
            break
        fi
    done
fi

LEVELS=("${@:-}")
if [ $# -eq 0 ]; then
    LEVELS=(level.0 level.1 level.2 level.3)
fi

PASS=0; FAIL=0; CRASH=0; ABORT=0; KNOWN=0; SKIP=0
FAILED_LIST=()

# tahoe exits 0 even when a run gives up (issue #74): after an unrecovered solver failure it prints
# "Exiting time sequence" (load-step cutting exhausted, run ends early), and after an exception it
# prints "exit on exception". Such runs are reported as ABORT and fail the job like a crash.
# Decks that are known to abort, each tracked by an issue, are listed (path relative to
# benchmark_XML, then "# #<issue> reason") in benchmark_XML/known_aborts.txt. They are reported as
# KNOWN and do not fail the job; any other abort does.
KNOWN_ABORTS_FILE="${KNOWN_ABORTS_FILE:-$BENCH_ROOT/known_aborts.txt}"
is_known_abort() {  # full path -> 0 if listed
    [ -f "$KNOWN_ABORTS_FILE" ] || return 1
    local rel="${1#"$BENCH_ROOT"/}"
    grep -v '^[[:space:]]*#' "$KNOWN_ABORTS_FILE" | awk '{print $1}' | grep -qxF "$rel"
}

abort_reason() {  # log -> prints a one-line reason, or nothing if the run completed
    local log="$1"
    if grep -q "exit on exception" "$log"; then
        grep -m1 -A1 "exit on exception" "$log" | tr '\n' ' ' | tr -s ' ' | cut -c1-120
    elif grep -q "Exiting time sequence" "$log"; then
        local step; step=$(grep -E "^ *Step: " "$log" | tail -1 | tr -s ' ' | sed 's/^ //')
        echo "time sequence ended early (load-step cutting exhausted) at ${step:-unknown step}"
    fi
}

# Process a single XML benchmark file (serial)
run_one() {
    local dir="$1"   # directory containing the xml
    local xml="$2"   # xml filename (basename only)
    local full="$dir/$xml"

    [ -f "$full" ] || { ((SKIP++)); return; }

    # Run tahoe
    local rc=0 log
    log=$(mktemp)
    (cd "$dir" && timeout 120 "$TAHOE" -f "$xml" > "$log" 2>&1) || rc=$?
    if [ $rc -eq 124 ]; then
        echo "  TIMEOUT  $full"
        ((CRASH++)); FAILED_LIST+=("TIMEOUT: $full"); rm -f "$log"; return
    elif [ $rc -ne 0 ]; then
        echo "  CRASH    $full  (exit $rc)"
        ((CRASH++)); FAILED_LIST+=("CRASH:   $full"); rm -f "$log"; return
    fi
    local why; why=$(abort_reason "$log"); rm -f "$log"
    if [ -n "$why" ]; then
        if is_known_abort "$full"; then
            echo "  KNOWN    $full  (listed in known_aborts.txt: $why)"
            ((KNOWN++)); return
        fi
        echo "  ABORT    $full  ($why)"
        ((ABORT++)); FAILED_LIST+=("ABORT:   $full"); return
    fi

    compare_one "$dir" "$xml" "$full"
}

# Compare a single XML result against benchmark reference
compare_one() {
    local dir="$1"
    local xml="$2"
    local full="$3"

    # Check if there are any reference files to compare against
    local stem="${xml%.xml}"
    if ! ls "$dir/benchmark/${stem}"* > /dev/null 2>&1; then
        # No benchmark reference — just count the run as pass (no comparison)
        ((PASS++)); return
    fi

    # Run comparator
    local cmp_out rc2=0
    cmp_out=$(cd "$dir" && timeout 30 "$COMPARE" -f "$xml" 2>&1) || rc2=$?

    if echo "$cmp_out" | grep -q ": PASS"; then
        ((PASS++))
    elif echo "$cmp_out" | grep -q ": FAIL"; then
        local detail
        detail=$(echo "$cmp_out" | grep -E "FAIL|differ" | head -3 | tr '\n' '; ')
        echo "  FAIL     $full"
        echo "           $detail"
        ((FAIL++)); FAILED_LIST+=("FAIL:    $full")
    else
        # comparator had an error
        echo "  ERROR    $full  (compare rc=$rc2)"
        ((FAIL++)); FAILED_LIST+=("ERROR:   $full")
    fi
}

# Recursively process a run.batch file (serial)
process_batch() {
    local dir="$1"
    local batch="$dir/run.batch"
    [ -f "$batch" ] || return 0

    while IFS= read -r line; do
        line="${line%%#*}"          # strip inline comments
        line="${line#"${line%%[![:space:]]*}"}"  # ltrim
        line="${line%"${line##*[![:space:]]}"}"  # rtrim
        [ -z "$line" ] && continue
        [ "${line:0:1}" = "@" ] && continue   # batch job-char line

        if [[ "$line" == *.xml ]]; then
            run_one "$dir" "$line"
        elif [[ "$line" == *run.batch ]]; then
            local subdir="$dir/${line%/run.batch}"
            process_batch "$subdir"
        fi
    done < "$batch"
}

# Collect all XML files referenced in a level.3 batch hierarchy
# Sets global array LEVEL3_XMLS as pairs "dir xml"
collect_level3_xmls() {
    local dir="$1"
    local batch="$dir/run.batch"
    [ -f "$batch" ] || return 0

    local skip_options=1
    while IFS= read -r line; do
        line="${line%%#*}"
        line="${line#"${line%%[![:space:]]*}"}"
        line="${line%"${line##*[![:space:]]}"}"
        [ -z "$line" ] && continue
        [ "${line:0:1}" = "@" ] && continue
        [ "${line:0:1}" = "-" ] && continue   # command-line option; skip

        if [[ "$line" == *.xml ]]; then
            LEVEL3_XMLS+=("$dir|$line")
        elif [[ "$line" == *run.batch ]]; then
            local subdir="$dir/${line%/run.batch}"
            collect_level3_xmls "$subdir"
        fi
    done < "$batch"
}

# Run level.3 MPI tests
run_level3() {
    local parallel_dir="$BENCH_ROOT/level.3/parallel"
    [ -d "$parallel_dir" ] || return

    if [ -z "$MPIRUN" ]; then
        echo "  SKIP (MPI-only)  $parallel_dir"
        echo "  System mpirun not found. Install openmpi: sudo apt-get install libopenmpi-dev"
        ((SKIP++))
        return
    fi

    echo "  Running MPI batch with $MPIRUN -np $MPI_NP ..."
    local rc=0 log
    log=$(mktemp)
    # shellcheck disable=SC2086
    (cd "$parallel_dir" && timeout 300 "$MPIRUN" --allow-run-as-root $MPIFLAGS \
        -np "$MPI_NP" "$TAHOE" -f run.batch > "$log" 2>&1) || rc=$?
    local why; why=$(abort_reason "$log"); rm -f "$log"
    if [ $rc -eq 0 ] && [ -n "$why" ]; then
        echo "  ABORT    level.3 MPI batch  ($why)"
        ((ABORT++)); FAILED_LIST+=("ABORT:   level.3/parallel/run.batch")
    fi
    if [ $rc -eq 124 ]; then
        echo "  TIMEOUT  level.3 MPI batch (300s)"
        ((CRASH++)); FAILED_LIST+=("TIMEOUT: level.3/parallel/run.batch"); return
    elif [ $rc -ne 0 ]; then
        echo "  CRASH    level.3 MPI batch  (exit $rc)"
        ((CRASH++)); FAILED_LIST+=("CRASH:   level.3/parallel/run.batch"); return
    fi

    # Collect and compare all XML files that were run
    LEVEL3_XMLS=()
    collect_level3_xmls "$parallel_dir"

    for entry in "${LEVEL3_XMLS[@]}"; do
        local dir="${entry%%|*}"
        local xml="${entry##*|}"
        local full="$dir/$xml"
        [ -f "$dir/$xml" ] || continue
        compare_one "$dir" "$xml" "$full"
    done
}

for level in "${LEVELS[@]}"; do
    level_dir="$BENCH_ROOT/$level"
    [ -d "$level_dir" ] || { echo "Skipping missing level: $level"; continue; }
    echo "========================================================"
    echo "  $level"
    echo "========================================================"

    if [ "$level" = "level.3" ]; then
        run_level3
    elif [ -f "$level_dir/run.batch" ]; then
        process_batch "$level_dir"
    else
        # No top-level run.batch — check subdirectories
        found=0
        for subdir in "$level_dir"/*/; do
            [ -f "$subdir/run.batch" ] && { process_batch "$subdir"; found=1; }
        done
        [ $found -eq 0 ] && echo "  (no run.batch found)"
    fi

    echo ""
done

echo "========================================================"
echo "  SUMMARY"
echo "========================================================"
printf "  %-8s %d\n" PASS    $PASS
printf "  %-8s %d\n" FAIL    $FAIL
printf "  %-8s %d\n" CRASH   $CRASH
printf "  %-8s %d\n" ABORT   $ABORT
printf "  %-8s %d   (known aborts, see known_aborts.txt)\n" KNOWN   $KNOWN
printf "  %-8s %d\n" SKIP    $SKIP
echo ""
if [ ${#FAILED_LIST[@]} -gt 0 ]; then
    echo "  FAILURES / ERRORS:"
    for f in "${FAILED_LIST[@]}"; do
        echo "    $f"
    done
fi

# Exit non-zero if any run crashed (non-zero exit or timeout) or aborted (tahoe gave up but exited 0,
# issue #74). Reference-comparison failures (FAIL) are tolerated — many are known issues with
# optional modules or ExodusII format mismatches.
if [ $CRASH -gt 0 ] || [ $ABORT -gt 0 ]; then
    echo ""
    echo "  ** $CRASH crash(es), $ABORT aborted run(s) detected — CI FAIL **"
    exit 1
fi
