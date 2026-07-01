#!/usr/bin/env python3
"""Elastocapillary sweep for mechanical compression: one run per elastocapillary
number gbar.  Two-stage loading (ramp surface tension to target over [0,n], then
linear displacement ramp [n,m] with small increment) on a long film so several
wrinkles fit (wavelength measurable).  BCs: bottom roller, both ends pushed inward
equally (symmetric compression -> no loaded-edge bulge), top free.  Explicit
central-difference, mass scaling + light damping.
Usage: python3 sweep.py <gbar> <eps_max> [tag]
"""
import os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
TAHOE = "/home/samanseifi/codes/tahoe/build/bin/tahoe"
MESH = os.environ.get("MESH", "film_L80.geom")
L = float(os.environ.get("L", 80.0))
H = float(os.environ.get("H", 4.0))
Nx = int(os.environ.get("NX", 160))
DT = float(os.environ.get("DT", 0.005))
DENS = float(os.environ.get("DENS", 5.0))
DAMP = float(os.environ.get("DAMP", 0.1))
# quasi-static strain rate scaled so the ramp = ~4 wave-crossings regardless of L
RATE = float(os.environ.get("RATE", 0.0007 * 80.0 / L))
N = float(os.environ.get("N", 50.0))     # surface-tension ramp time (t_0); displacement held to t=N
NFRAMES = int(os.environ.get("NFRAMES", 150))   # output frames (dense -> resolves sharp onset)


def deck(gbar, eps_max, ramp, vhalf, nsteps, out_inc):
    gamma = gbar * H
    M = N + ramp
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<tahoe author="sam" compute_IC="false" geometry_file="{MESH}"
    output_format="ExodusII" title="elastocapillary sweep gbar={gbar}"
    xmlns:x0="http://www.w3.org/2001/XMLSchema">
  <time num_steps="{nsteps}" output_inc="{out_inc}" time_step="{DT}">
    <schedule_function><piecewise_linear>
      <OrderedPair x="0" y="0"/><OrderedPair x="{N}" y="0"/>
      <OrderedPair x="{M}" y="{ramp}"/><OrderedPair x="10000000" y="{ramp}"/>
    </piecewise_linear></schedule_function>
  </time>
  <nodes>
    <field field_name="displacement" integrator="central_difference" solution_group="1">
      <dof_labels><String value="D_X"/><String value="D_Y"/></dof_labels>
      <kinematic_BC dof="2" node_ID="1" type="fixed"/>                                 <!-- bottom roller u_y=0 -->
      <kinematic_BC dof="1" node_ID="3" schedule="1" type="u" value="{vhalf:.6f}"/>     <!-- left pushed inward (+x) -->
      <kinematic_BC dof="1" node_ID="4" schedule="1" type="u" value="{-vhalf:.6f}"/>    <!-- right pushed inward (-x) -->
    </field>
  </nodes>
  <element_list>
    <updated_lagrangian_Q1P0_surface field_name="displacement" mass_type="lumped_mass" epsilon="1.0">
      <quadrilateral/><solid_element_nodal_output displacements="1" stress="1"/>
      <large_strain_element_block><block_ID_list><String value="1"/></block_ID_list>
        <large_strain_material_2D><RG_split_general constraint_2D="plane_strain" density="{DENS}">
          <rg_eq_potential><neo-hookean kappa="1000.67" mu="1.0"/></rg_eq_potential>
        </RG_split_general></large_strain_material_2D></large_strain_element_block>
      <surface_tension side_set_ID="1" gamma="{gamma}" t_0="{N}"/>
    </updated_lagrangian_Q1P0_surface>
  </element_list>
  <linear_solver><diagonal_matrix/></linear_solver>
  <solver_phases max_loops="1"><solver_phase iterations="1" pass_iterations="0" solver="1"/></solver_phases>
</tahoe>
"""


def run(gbar, eps_max, tag=None):
    ramp = eps_max / RATE
    vhalf = eps_max * L / (2.0 * ramp)            # each end pushed inward by eps*L/2 (symmetric)
    nsteps = int((N + ramp + 100) / DT)
    out_inc = max(1, nsteps // NFRAMES)
    name = tag or f"sweep_g{int(round(gbar*100)):04d}"
    open(os.path.join(HERE, name + ".xml"), "w").write(
        deck(gbar, eps_max, ramp, vhalf, nsteps, out_inc))
    env = dict(os.environ, TAHOE_DAMP=str(DAMP))
    t0 = time.time()
    rc = subprocess.run([TAHOE, "-f", name + ".xml"], cwd=HERE, env=env,
                        stdout=open(os.path.join(HERE, name + ".log"), "w"),
                        stderr=subprocess.STDOUT).returncode
    print(f"{name}: gbar={gbar} eps_max={eps_max} rc={rc} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    gbar = float(sys.argv[1]); eps_max = float(sys.argv[2])
    tag = sys.argv[3] if len(sys.argv) > 3 else None
    run(gbar, eps_max, tag)
