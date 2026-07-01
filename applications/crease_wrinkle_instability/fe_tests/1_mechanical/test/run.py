#!/usr/bin/env python3
"""Mechanical compression of a fine film -> ParaView (ExodusII).
Two-stage loading with the (0,0)->(n,0)->(m,m-n) displacement schedule:
  stage 1  t in [0, n] : surface tension ramps to target (t_0 = n); displacement held at 0
  stage 2  t in [n, m] : displacement ramps linearly (schedule slope = 1).
The schedule value at t=m is (m-n), so with kinematic_BC value = -eps*L/(m-n) the
total push is -eps*L and the per-step displacement increment value*1*dt is tiny.
BCs: bottom roller (u_y=0), left fixed in x, right pushed in, top free.
Explicit central-difference, mass scaling (density), light damping (TAHOE_DAMP).
Open compress.io0.exo in ParaView.
"""
import os, subprocess, time

HERE = os.path.dirname(os.path.abspath(__file__))
TAHOE = "/home/samanseifi/codes/tahoe/build/bin/tahoe"
L, H = 40.0, 4.0
EPS = float(os.environ.get("EPS", 0.60))      # final compressive strain
GAMMA = float(os.environ.get("GAMMA", 4.0))   # target surface tension (gbar = GAMMA/(mu*H))
N = float(os.environ.get("N", 50.0))          # surface-tension ramp time (= t_0); disp held to t=N
RAMP = float(os.environ.get("RAMP", 900.0))   # m-n: displacement ramp duration (larger = slower, smaller increment)
DT = 0.005
DENS = float(os.environ.get("DENS", 5.0))     # mild mass scaling (DT stability)
DAMP = float(os.environ.get("DAMP", 0.1))     # light mass-proportional damping
M = N + RAMP
push = EPS * L
value = -push / RAMP                           # small per-step displacement increment (= value*dt)
NSTEPS = int((M + 100) / DT)
NFRAMES = int(os.environ.get("NFRAMES", 100))  # target number of ParaView frames
OUTPUT_INC = max(1, NSTEPS // NFRAMES)

DECK = f"""<?xml version="1.0" encoding="UTF-8"?>
<tahoe author="sam" compute_IC="false" geometry_file="film.geom"
    output_format="ExodusII" title="film compression (ST ramp then displacement)"
    xmlns:x0="http://www.w3.org/2001/XMLSchema">
  <time num_steps="{NSTEPS}" output_inc="{OUTPUT_INC}" time_step="{DT}">
    <schedule_function><piecewise_linear>
      <OrderedPair x="0" y="0"/>
      <OrderedPair x="{N}" y="0"/>
      <OrderedPair x="{M}" y="{RAMP}"/>
      <OrderedPair x="10000000" y="{RAMP}"/>
    </piecewise_linear></schedule_function>
  </time>
  <nodes>
    <field field_name="displacement" integrator="central_difference" solution_group="1">
      <dof_labels><String value="D_X"/><String value="D_Y"/></dof_labels>
      <kinematic_BC dof="2" node_ID="1" type="fixed"/>                                 <!-- bottom roller u_y=0 -->
      <kinematic_BC dof="1" node_ID="3" type="fixed"/>                                 <!-- left fixed in x -->
      <kinematic_BC dof="1" node_ID="4" schedule="1" type="u" value="{value:.6f}"/>    <!-- right pushed (small increment) -->
    </field>
  </nodes>
  <element_list>
    <updated_lagrangian_Q1P0_surface field_name="displacement" mass_type="lumped_mass" epsilon="1.0">
      <quadrilateral/><solid_element_nodal_output displacements="1" stress="1"/>
      <large_strain_element_block><block_ID_list><String value="1"/></block_ID_list>
        <large_strain_material_2D><RG_split_general constraint_2D="plane_strain" density="{DENS}">
          <rg_eq_potential><neo-hookean kappa="1000.67" mu="1.0"/></rg_eq_potential>
        </RG_split_general></large_strain_material_2D></large_strain_element_block>
      <surface_tension side_set_ID="1" gamma="{GAMMA}" t_0="{N}"/>     <!-- stage 1: ramp to target over [0,n] -->
    </updated_lagrangian_Q1P0_surface>
  </element_list>
  <linear_solver><diagonal_matrix/></linear_solver>
  <solver_phases max_loops="1"><solver_phase iterations="1" pass_iterations="0" solver="1"/></solver_phases>
</tahoe>
"""


def main():
    open(os.path.join(HERE, "compress.xml"), "w").write(DECK)
    env = dict(os.environ, TAHOE_DAMP=str(DAMP))
    print(f"stage1: gamma->{GAMMA} (gbar={GAMMA/(1*H):.2f}) over t=[0,{N:.0f}]; "
          f"stage2: eps->{EPS} over t=[{N:.0f},{M:.0f}]  (disp value={value:.5f}/unit, "
          f"increment {abs(value)*DT:.2e}/step), {NSTEPS} steps, damp={DAMP}", flush=True)
    t0 = time.time()
    rc = subprocess.run([TAHOE, "-f", "compress.xml"], cwd=HERE, env=env,
                        stdout=open(os.path.join(HERE, "compress.log"), "w"),
                        stderr=subprocess.STDOUT).returncode
    print(f"done: rc={rc} ({time.time()-t0:.0f}s) -> open compress.io0.exo in ParaView", flush=True)


if __name__ == "__main__":
    main()
