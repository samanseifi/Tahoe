#!/usr/bin/env python3
"""Electromechanical instability sweep, UNSTRAINED film (lambda_pre = 1), UNDAMPED.
Staggered two-field scheme: electric_scalar_potential (static, MUMPS) + displacement
(explicit central difference).  A flat film is laterally confined (u_x=0 on both
sides), grounded on the bottom, and a voltage is ramped on the compliant top
electrode until the surface wrinkles under the Maxwell traction.  Surface tension
gamma = gbar*H_f on the top.  No artificial damping (like the mechanical campaign)
so the wrinkle develops in the bulk.
  Stage 1 (t in [0,N]):   ramp surface tension (t_0=N), voltage held 0.
  Stage 2 (t in [N,M]):   ramp voltage 0 -> V_max.
Measure: onset (top-surface amplitude), critical nominal field Etilde_c = Phi_c/H_f
(normalized Etilde_c*sqrt(eps/mu) with eps=mu=1), and wavelength.
Usage: python3 em_sweep.py <gbar> [tag]
"""
import os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
TAHOE = "/home/samanseifi/codes/tahoe/build/bin/tahoe"
MESH = "film_L120.geom"
L, H = 120.0, 4.0
EPS_PERM = 1.0                     # dielectric permittivity
MU = 1.0
gbar = float(sys.argv[1])
tag = sys.argv[2] if len(sys.argv) > 2 else f"em_g{int(round(gbar*100)):04d}"
gamma = gbar * H
VMAX = float(os.environ.get("VMAX", 18.0))     # top-electrode voltage; Etilde=V/H up to 4.5
N = float(os.environ.get("N", 50.0))           # surface-tension ramp (t_0); voltage held to t=N
RAMP = float(os.environ.get("RAMP", 1000.0))   # voltage ramp duration
M = N + RAMP
DT = float(os.environ.get("DT", 0.005))
DENS = float(os.environ.get("DENS", 5.0))
DAMP = float(os.environ.get("DAMP", 0.0))      # UNDAMPED by default
NSTEPS = int((M + 100) / DT)
NFRAMES = int(os.environ.get("NFRAMES", 300))
OUT_INC = max(1, NSTEPS // NFRAMES)

DECK = f"""<?xml version="1.0" encoding="UTF-8"?>
<tahoe author="sam" compute_IC="false" geometry_file="{MESH}" output_format="ExodusII"
    title="EM unstrained gbar={gbar}" xmlns:x0="http://www.w3.org/2001/XMLSchema">
  <time num_steps="{NSTEPS}" output_inc="{OUT_INC}" time_step="{DT}">
    <schedule_function><piecewise_linear>
      <OrderedPair x="0" y="0"/><OrderedPair x="{N}" y="0"/>
      <OrderedPair x="{M}" y="1"/><OrderedPair x="10000000" y="1"/>
    </piecewise_linear></schedule_function>
  </time>
  <nodes>
    <field field_name="electric_scalar_potential" solution_group="1" integrator="static">
      <dof_labels><String value="Psi"/></dof_labels>
      <kinematic_BC node_ID="1" dof="1" type="fixed" value="0.0"/>              <!-- bottom ground -->
      <kinematic_BC node_ID="2" dof="1" type="u" schedule="1" value="{VMAX}"/>  <!-- top electrode -->
    </field>
    <field field_name="displacement" solution_group="2" integrator="central_difference">
      <dof_labels><String value="D_X"/><String value="D_Y"/></dof_labels>
      <kinematic_BC node_ID="1" dof="2" type="fixed" value="0.0"/>              <!-- bottom roller u_y=0 -->
      <kinematic_BC node_ID="3" dof="1" type="fixed" value="0.0"/>              <!-- left  u_x=0 (confined) -->
      <kinematic_BC node_ID="4" dof="1" type="fixed" value="0.0"/>              <!-- right u_x=0 (confined) -->
    </field>
  </nodes>
  <element_list>
    <diffusion field_name="electric_scalar_potential">
      <quadrilateral num_ip="4"/>
      <diffusion_element_nodal_output coordinates="1" displacement="1"/>
      <diffusion_element_block><block_ID_list><String value="1"/></block_ID_list>
        <diffusion_material><linear_dielectric_material epsilon="{EPS_PERM}"/></diffusion_material>
      </diffusion_element_block>
    </diffusion>
    <updated_lagrangian_Q1P0_surface field_name="displacement" mass_type="lumped_mass" epsilon="1.0">
      <quadrilateral num_ip="4"/><solid_element_nodal_output displacements="1" stress="1"/>
      <large_strain_element_block><block_ID_list><String value="1"/></block_ID_list>
        <large_strain_material_2D><RG_split_general constraint_2D="plane_strain" density="{DENS}">
          <rg_eq_potential><neo-hookean kappa="1000.67" mu="{MU}"/></rg_eq_potential>
        </RG_split_general></large_strain_material_2D></large_strain_element_block>
      <surface_tension side_set_ID="1" gamma="{gamma}" t_0="{N}"/>
    </updated_lagrangian_Q1P0_surface>
  </element_list>
  <nonlinear_solver abs_tolerance="1.0e-10" rel_tolerance="1.0e-10" max_iterations="15" divergence_tolerance="1.0e8">
    <MUMPS_matrix message_level="silent" always_symmetric="false"/>
  </nonlinear_solver>
  <linear_solver><diagonal_matrix/></linear_solver>
  <solver_phases max_loops="1">
    <solver_phase solver="1" iterations="1" pass_iterations="0"/>
    <solver_phase solver="2" iterations="1" pass_iterations="0"/>
  </solver_phases>
</tahoe>
"""


def main():
    open(os.path.join(HERE, tag + ".xml"), "w").write(DECK)
    env = dict(os.environ, TAHOE_DAMP=str(DAMP))
    print(f"{tag}: gbar={gbar} gamma={gamma} Vmax={VMAX} (Etilde_max={VMAX/H:.2f}) "
          f"undamped={DAMP==0} {NSTEPS} steps", flush=True)
    t0 = time.time()
    rc = subprocess.run([TAHOE, "-f", tag + ".xml"], cwd=HERE, env=env,
                        stdout=open(os.path.join(HERE, tag + ".log"), "w"),
                        stderr=subprocess.STDOUT).returncode
    print(f"{tag}: rc={rc} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
