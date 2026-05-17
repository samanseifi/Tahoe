#!/bin/bash
# Benchmark: 3D surface tension, SPOOLES (1p) vs MUMPS MPI (1/2/4/6 cores)
TAHOE=/home/samanseifi/codes/tahoe/build-mpi-mumps/bin/tahoe
cd /home/samanseifi/codes/tahoe/benchmark_XML/level.4/surface_tension

run_and_time() {
    local label="$1"; local xml="$2"; shift 2; local cmd=("$@")
    rm -f "${xml%.xml}".out "${xml%.xml}".p*.out "${xml%.xml}".p*.log \
          "${xml%.xml}".io*.exo "${xml%.xml}".echo.xml "${xml%.xml}".valid.xml
    t0=$(date +%s%N)
    "${cmd[@]}" -f "$xml" > /dev/null 2>&1
    t1=$(date +%s%N)
    wall=$(echo "scale=2; ($t1 - $t0)/1000000000" | bc)
    # Timing from tahoe's rank-0 .out file
    outf="${xml%.xml}.out"
    [ ! -f "$outf" ] && outf="${xml%.xml}.p0.out"
    cons=$(grep "Construction:" "$outf" 2>/dev/null | awk '{print $2}')
    solv=$(grep "Solution:"     "$outf" 2>/dev/null | awk '{print $2}')
    printf "  %-26s  wall=%6.2fs  construct=%-8s  solve=%-8s\n" \
        "$label" "$wall" "${cons:-N/A}" "${solv:-N/A}"
}

echo "================================================================"
echo " 3D Surface Tension Benchmark"
echo " Mesh: block_3D — 1331 nodes, 1000 hex, 3993 DOFs"
echo " Steps: 100 (no output, gamma=5.1, t_0=50)"
echo "================================================================"
echo ""
echo "SPOOLES (serial):"
run_and_time "SPOOLES 1-process" test_3D_bench_spooles.xml \
    mpirun -np 1 "$TAHOE"

echo ""
echo "MUMPS MPI:"
for np in 1 2 4 6; do
    run_and_time "MUMPS ${np}-process(es)" test_3D_bench_mumps.xml \
        mpirun -np "$np" "$TAHOE"
done

echo ""
echo "================================================================"
echo " Duplicate-error-print check (mpirun -np 2, count 'Absolute error' lines)"
echo "================================================================"
mpirun -np 2 "$TAHOE" -f test_3D_bench_mumps.xml 2>&1 | \
    grep -c "Absolute error" | xargs -I{} echo "  Lines printed: {} (expect 1 per step, ~100 total — NOT doubled to ~200)"
