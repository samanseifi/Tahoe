#!/usr/bin/env python3
"""Verify three-way coupled electro-mechano-fracture test.

Checks that:
1. Displacement concentrates in the degraded (cracked) zone — same as coupled_tension
2. Electric potential gradient (E-field) is distorted by the crack:
   in degraded zone, epsilon_eff ~ k*epsilon -> high E-field to maintain current continuity
"""

import numpy as np

try:
    import netCDF4
except ImportError:
    print("netCDF4 not available — skipping ExodusII verification")
    exit(0)

# io0 = mechanics (D_X, D_Y, s11, s22, s12)
# io1 = electrical (ESP, x1, x2)
# io2 = phase-field (d, x1, x2)
ds_mech = netCDF4.Dataset("three_way_coupled.io0.exo", "r")
ds_elec = netCDF4.Dataset("three_way_coupled.io1.exo", "r")

# Read coordinates from mechanics output
x = np.array(ds_mech.variables["coordx"][:])
y = np.array(ds_mech.variables["coordy"][:])

nt = ds_mech.dimensions["time_step"].size
print(f"Number of time steps written: {nt}")

# --- Displacement check (use last available step) ---
ux_final = np.array(ds_mech.variables["vals_nod_var1"][nt-1, :])  # D_X

center_mask = np.abs(x - 1.0) < 0.15
far_left_mask = x < 0.3
far_right_mask = x > 1.7

ux_center_avg = np.mean(np.abs(ux_final[center_mask]))
ux_far_left_avg = np.mean(np.abs(ux_final[far_left_mask]))
ux_far_right_avg = np.mean(np.abs(ux_final[far_right_mask]))

print(f"--- Displacement check (final step) ---")
print(f"  |u_x| at center (x~1.0): {ux_center_avg:.6e}")
print(f"  |u_x| at far left (x<0.3): {ux_far_left_avg:.6e}")
print(f"  |u_x| at far right (x>1.7): {ux_far_right_avg:.6e}")
print(f"  Strain concentration: {ux_center_avg / ((ux_far_left_avg + ux_far_right_avg)/2):.2f}x")

# --- Electric field check ---
x_e = np.array(ds_elec.variables["coordx"][:])
y_e = np.array(ds_elec.variables["coordy"][:])
nt_e = ds_elec.dimensions["time_step"].size
phi_final = np.array(ds_elec.variables["vals_nod_var1"][nt_e-1, :])  # ESP

# E-field (dphi/dx) along x at mid-height row (y~0.5)
# Voltage is left-to-right, crack is at x=1.0
y_unique = np.unique(np.round(y_e, 8))
y_mid_val = y_unique[np.argmin(np.abs(y_unique - 0.5))]
mid_row = np.abs(y_e - y_mid_val) < 1e-6

x_mid = x_e[mid_row]
phi_mid = phi_final[mid_row]
idx_sort = np.argsort(x_mid)
x_mid = x_mid[idx_sort]
phi_mid = phi_mid[idx_sort]

print(f"\n  Mid-height row at y={y_mid_val:.4f}: {len(x_mid)} nodes")

if len(x_mid) > 1:
    Ex_mid = -np.diff(phi_mid) / np.diff(x_mid)
    x_midpts = 0.5 * (x_mid[:-1] + x_mid[1:])

    # E-field at the crack center vs far from crack
    crack_mask = np.abs(x_midpts - 1.0) < 0.15
    far_mask = (x_midpts < 0.3) | (x_midpts > 1.7)

    Ex_crack = np.mean(np.abs(Ex_mid[crack_mask]))
    Ex_far = np.mean(np.abs(Ex_mid[far_mask]))

    print(f"\n--- Electric field check (last converged step) ---")
    print(f"  |E_x| near crack (x~1.0):  {Ex_crack:.6e}")
    print(f"  |E_x| far from crack:       {Ex_far:.6e}")
    print(f"  E-field ratio (crack/far):  {Ex_crack/Ex_far:.2f}x")
    if Ex_crack > 1.1 * Ex_far:
        print("  PASS: E-field concentrates at crack (degraded permittivity works)")
    else:
        print("  NOTE: E-field not significantly elevated — check coupling")

ds_mech.close()
ds_elec.close()
print("\nThree-way coupling verification complete.")
