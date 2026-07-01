#!/usr/bin/env python3
"""Periodic-BC compression test: no ends -> uniform bulk wrinkle (no Saint-Venant
edge effect).  Left face (NS3, leader) <-> right face (NS4, follower) are tied by
MappedPeriodicT (<mapped_nodes>): u_follower = u_leader + s(t) F_perturb (X_f-X_l),
with F_perturb=[[-eps,0],[0,0]] so u_x(right)=u_x(left)-s eps L (compression) and
u_y(right)=u_y(left) (periodic wrinkle).  NS5 (corner at origin) is affinely mapped
-> pinned at u=0 (rigid-body pin, and satisfies the required mapped_node_ID_list).
Bottom NS1 roller (u_y=0), top NS2 free + surface tension.
Usage: python3 run_periodic.py <gbar> <eps_max> [L] [tag]
"""
import os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
TAHOE = "/home/samanseifi/codes/tahoe/build/bin/tahoe"
gbar = float(sys.argv[1]); EPS = float(sys.argv[2])
L = float(sys.argv[3]) if len(sys.argv) > 3 else 72.0
tag = sys.argv[4] if len(sys.argv) > 4 else f"per_g{int(round(gbar*100)):04d}"
H = 4.0; MESH = f"film_L{int(L)}.geom"
gamma = gbar * H
N = 50.0; RAMP = float(os.environ.get("RAMP", 1000.0)); M = N + RAMP
DT = 0.005; DENS = 5.0; DAMP = float(os.environ.get("DAMP", 0.1))
NSTEPS = int((M + 100) / DT); NFRAMES = int(os.environ.get("NFRAMES", 150))
OUT_INC = max(1, NSTEPS // NFRAMES)

DECK = f"""<?xml version="1.0" encoding="UTF-8"?>
<tahoe author="sam" compute_IC="false" geometry_file="{MESH}"
    output_format="ExodusII" title="periodic compression gbar={gbar}"
    xmlns:x0="http://www.w3.org/2001/XMLSchema">
  <time num_steps="{NSTEPS}" output_inc="{OUT_INC}" time_step="{DT}">
    <schedule_function><piecewise_linear>
      <OrderedPair x="0" y="0"/><OrderedPair x="{N}" y="0"/>
      <OrderedPair x="{M}" y="1"/><OrderedPair x="10000000" y="1"/>
    </piecewise_linear></schedule_function>
  </time>
  <nodes>
    <field field_name="displacement" integrator="central_difference" solution_group="1">
      <dof_labels><String value="D_X"/><String value="D_Y"/></dof_labels>
      <kinematic_BC dof="2" node_ID="1" type="fixed"/>                <!-- bottom roller u_y=0 -->
      <mapped_nodes schedule="1">
        <Matrix_2x2 A_1_1="{-EPS}" A_1_2="0.0" A_2_1="0.0" A_2_2="0.0"/>
        <mapped_node_ID_list><String value="5"/></mapped_node_ID_list>
        <leader_node_ID_list><String value="3"/></leader_node_ID_list>
        <follower_node_ID_list><String value="4"/></follower_node_ID_list>
      </mapped_nodes>
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

open(os.path.join(HERE, tag + ".xml"), "w").write(DECK)
env = dict(os.environ, TAHOE_DAMP=str(DAMP))
print(f"{tag}: gbar={gbar} eps={EPS} L={L} periodic, {NSTEPS} steps, damp={DAMP}", flush=True)
t0 = time.time()
rc = subprocess.run([TAHOE, "-f", tag + ".xml"], cwd=HERE, env=env,
                    stdout=open(os.path.join(HERE, tag + ".log"), "w"),
                    stderr=subprocess.STDOUT).returncode
print(f"{tag}: rc={rc} ({time.time()-t0:.0f}s)", flush=True)
