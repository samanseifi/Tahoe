#!/usr/bin/env python3
"""Verify the phase-field screened Poisson solution against the analytical result.

Analytical solution: d(x) = exp(-|x|/ell)

This script reads the Tahoe output and compares nodal d values to the exact solution.
Usage:
    python3 verify_solution.py <output_file.exo>

If no ExodusII file is provided, it reads the .run output file for nodal values.
"""

import numpy as np
import sys

def analytical_solution(x, ell=1.0):
    """Exact 1D solution: d(x) = exp(-|x|/ell)"""
    return np.exp(-np.abs(x) / ell)

def read_geom_file(filename="bar_phase_field.geom"):
    """Read node coordinates from Tahoe geometry file."""
    coords = {}
    with open(filename) as f:
        lines = f.readlines()

    # Find *nodes section
    in_nodes = False
    num_nodes = 0
    nsd = 0
    count = 0
    for line in lines:
        line = line.strip()
        if line == "*nodes":
            in_nodes = True
            continue
        if in_nodes:
            if num_nodes == 0:
                num_nodes = int(line.split('#')[0].strip())
                continue
            if nsd == 0:
                nsd = int(line.split('#')[0].strip())
                continue
            if line == '' or line.startswith('#'):
                continue
            parts = line.split()
            nid = int(parts[0])
            x = float(parts[1])
            y = float(parts[2]) if nsd > 1 else 0.0
            coords[nid] = (x, y)
            count += 1
            if count >= num_nodes:
                break
    return coords

def main():
    # Read geometry for node coordinates
    coords = read_geom_file()

    # Parameters
    ell = 1.0
    L = 5.0

    # For now, just print what the expected solution should be
    # at the bottom row of nodes (y=0)
    bottom_nodes = {nid: c for nid, c in coords.items() if abs(c[1]) < 1e-10}

    # Sort by x-coordinate
    sorted_nodes = sorted(bottom_nodes.items(), key=lambda item: item[1][0])

    print("Phase-field screened Poisson verification")
    print("=" * 60)
    print(f"Parameters: Gc=1.0, ell={ell}, L={L}")
    print(f"Analytical: d(x) = exp(-|x|/{ell})")
    print(f"Number of bottom-row nodes: {len(sorted_nodes)}")
    print()
    print(f"{'Node':>6s}  {'x':>10s}  {'d_exact':>12s}")
    print("-" * 35)

    x_vals = []
    d_exact = []
    for nid, (x, y) in sorted_nodes:
        d = analytical_solution(x, ell)
        x_vals.append(x)
        d_exact.append(d)
        # Print a subset
        if abs(x) <= 3.0 and abs(x * 10) % 5 < 1.0:
            print(f"{nid:6d}  {x:10.4f}  {d:12.6e}")

    print()
    print("Key values:")
    print(f"  d(0)   = {analytical_solution(0, ell):.6f}  (should be 1.0, imposed BC)")
    print(f"  d(ell) = {analytical_solution(ell, ell):.6f}  (should be ~0.3679)")
    print(f"  d(2*ell) = {analytical_solution(2*ell, ell):.6f}  (should be ~0.1353)")
    print(f"  d(L)   = {analytical_solution(L, ell):.6e}  (should be ~{np.exp(-L/ell):.4e})")

    # Try to read ExodusII output if available
    try:
        import exodus
        if len(sys.argv) > 1:
            exo_file = sys.argv[1]
        else:
            # Try to find the output file
            import glob
            exo_files = glob.glob("screened_poisson.exo") + glob.glob("screened_poisson*.e")
            if not exo_files:
                print("\nNo ExodusII output file found. Run tahoe first.")
                return
            exo_file = exo_files[0]

        print(f"\nReading ExodusII output: {exo_file}")
        e = exodus.exodus(exo_file, mode='r')
        # Read nodal variable 'd' at last time step
        num_times = e.num_times()
        d_names = e.get_node_variable_names()
        print(f"  Time steps: {num_times}")
        print(f"  Nodal variables: {d_names}")

        if 'd' in d_names:
            d_vals = e.get_node_variable_values('d', num_times)
            # Compare
            errors = []
            for i, (nid, (x, y)) in enumerate(sorted_nodes):
                d_num = d_vals[nid - 1]
                d_ex = analytical_solution(x, ell)
                errors.append(abs(d_num - d_ex))

            max_err = max(errors)
            l2_err = np.sqrt(np.mean(np.array(errors)**2))
            print(f"\n  Max error: {max_err:.6e}")
            print(f"  L2 error:  {l2_err:.6e}")
            if max_err < 0.01:
                print("  PASS: Solution matches analytical result.")
            else:
                print("  FAIL: Solution does not match analytical result.")

        e.close()
    except ImportError:
        print("\nNote: Install 'exodus' Python module to automatically verify output.")
        print("Otherwise, compare nodal 'd' values to d(x) = exp(-|x|/ell) manually.")

if __name__ == "__main__":
    main()
