#!/bin/bash
# 3D Benchmark: legacy vs vectorized explicit solver (Hex8 elements)

TAHOE=/home/samanseifi/codes/tahoe/build/bin/tahoe
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"

echo "============================================================"
echo "  3D Explicit Solver Benchmark: Legacy vs MVSIZ (Hex8)"
echo "  OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}"
echo "============================================================"
echo ""

for SIZE in small medium large; do
    LEGACY="legacy_hex_${SIZE}.xml"
    VECT="vectorized_hex_${SIZE}.xml"
    VECT1="vectorized_hex_${SIZE}_1ip.xml"

    echo "--- $SIZE mesh ---"

    # Legacy (8-IP)
    rm -f *.exo *.out *.echo.xml *.valid.xml SPOOLES.out
    echo -n "  Legacy 8-IP:       "
    T0=$(date +%s%N)
    $TAHOE -f "$LEGACY" > /dev/null 2>&1
    T1=$(date +%s%N)
    LEGACY_MS=$(( (T1 - T0) / 1000000 ))
    echo "${LEGACY_MS} ms"

    # Vectorized (8-IP)
    rm -f *.exo *.out *.echo.xml *.valid.xml SPOOLES.out
    echo -n "  Vectorized 8-IP:   "
    T0=$(date +%s%N)
    $TAHOE -f "$VECT" > /dev/null 2>&1
    T1=$(date +%s%N)
    VECT_MS=$(( (T1 - T0) / 1000000 ))
    echo "${VECT_MS} ms"

    # Vectorized (1-IP + hourglass)
    if [ -f "$VECT1" ]; then
        rm -f *.exo *.out *.echo.xml *.valid.xml SPOOLES.out
        echo -n "  Vectorized 1-IP+HG:"
        T0=$(date +%s%N)
        $TAHOE -f "$VECT1" > /dev/null 2>&1
        T1=$(date +%s%N)
        VECT1_MS=$(( (T1 - T0) / 1000000 ))
        echo " ${VECT1_MS} ms"
    fi

    # Speedups
    if [ $VECT_MS -gt 0 ]; then
        echo "  Speedup (8-IP): $(echo "scale=1; $LEGACY_MS / $VECT_MS" | bc)x"
    fi
    if [ -f "$VECT1" ] && [ $VECT1_MS -gt 0 ]; then
        echo "  Speedup (1-IP): $(echo "scale=1; $LEGACY_MS / $VECT1_MS" | bc)x"
    fi
    echo ""

    rm -f *.exo *.out *.echo.xml *.valid.xml SPOOLES.out
done
