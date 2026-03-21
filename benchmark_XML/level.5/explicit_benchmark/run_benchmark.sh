#!/bin/bash
# Benchmark: legacy vs vectorized explicit solver
# Runs both on small/medium/large meshes and reports wall-clock time.

TAHOE=/home/samanseifi/codes/tahoe/build/bin/tahoe
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"

echo "============================================"
echo "  Explicit Solver Benchmark: Legacy vs MVSIZ"
echo "============================================"
echo ""

for SIZE in small medium large; do
    LEGACY="legacy_${SIZE}.xml"
    VECT="vectorized_${SIZE}.xml"

    if [ ! -f "$LEGACY" ] || [ ! -f "$VECT" ]; then
        echo "Skipping $SIZE (XML files not found)"
        continue
    fi

    echo "--- $SIZE mesh ---"

    # Clean old output
    rm -f *.exo *.out *.echo.xml *.valid.xml SPOOLES.out

    # Legacy
    echo -n "  Legacy (updated_lagrangian): "
    T0=$(date +%s%N)
    $TAHOE -f "$LEGACY" > /dev/null 2>&1
    T1=$(date +%s%N)
    LEGACY_MS=$(( (T1 - T0) / 1000000 ))
    echo "${LEGACY_MS} ms"

    rm -f *.exo *.out *.echo.xml *.valid.xml SPOOLES.out

    # Vectorized
    echo -n "  Vectorized (explicit_solid): "
    T0=$(date +%s%N)
    $TAHOE -f "$VECT" > /dev/null 2>&1
    T1=$(date +%s%N)
    VECT_MS=$(( (T1 - T0) / 1000000 ))
    echo "${VECT_MS} ms"

    if [ $VECT_MS -gt 0 ]; then
        SPEEDUP=$(echo "scale=2; $LEGACY_MS / $VECT_MS" | bc)
        echo "  Speedup: ${SPEEDUP}x"
    fi
    echo ""

    rm -f *.exo *.out *.echo.xml *.valid.xml SPOOLES.out
done
