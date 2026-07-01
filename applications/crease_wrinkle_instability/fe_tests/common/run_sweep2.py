#!/usr/bin/env python3
"""
Extended electromechanical FEA sweep on a REFINED mesh.

Same two-stage workflow and identical voltage ramp as run_sweep.py (so extract.py
applies unchanged), but parameterized by:
  * MESH      : geometry file (default refined 160x16)
  * the gamma sweep at lambda=1:   GBARS  (elastocapillary up to 20)
  * the pre-strain sweep at a fixed gbar:  LAMS  (compression and stretch)

Ramp (unchanged):  gamma ramped by t=40; nominal field Etilde = Phi/H ramps
Etilde = 10*clip((t-40)/1200,0,1)  (PHI_MAX=40, H=4 -> Etilde_max=10).

Usage:
  python3 run_sweep2.py stage1            # static pre-strain restarts
  python3 run_sweep2.py stage2-decks      # write all stage2 decks (no run)
  python3 run_sweep2.py run <deck.xml>    # run one deck
  STEPS=4000 python3 run_sweep2.py time   # quick timing run (lam=1,gbar=1)
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TAHOE = "/home/samanseifi/codes/tahoe/build/bin/tahoe"
MESH = os.environ.get("MESH", "../../meshes/bar_2D_ref160x16.geom")
H = 4.0
PHI_MAX = 40.0
DT = float(os.environ.get("DT", 0.02))   # refined mesh (dy=0.25) needs DT<=~0.008 (CFL, kappa=1000)
T_GAMMA = 40.0
T_V0 = 40.0
T_V1 = 1240.0
NSTEPS = int(os.environ.get("STEPS", int(T_V1 / DT) + 2000))

# half-strip grip is +-(1-lam)*40 for the L=80 strip (node sets 3,4 are the ends)
GRIP = 40.0

# gamma sweep at lambda = 1
GBARS = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
# pre-strain sweep at fixed gbar
LAM_GBAR = 2.0
LAMS = [0.80, 0.85, 0.90, 0.95, 1.05, 1.10, 1.20, 1.30]


def lam_tag(lam):
    return f"{int(round(lam*100)):03d}"


def stage1_xml(lam):
    g3, g4 = (1.0 - lam) * GRIP, -(1.0 - lam) * GRIP
    # stage 1 is STATIC (implicit Newton): use a fixed step so 2000 steps reach t=40,
    # fully applying the grip schedule (full at t=39); independent of the explicit DT.
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<tahoe author="sam" compute_IC="false" geometry_file="{MESH}"
    output_format="ExodusII" restart_output_inc="2000" title="stage1 lam={lam}"
    xmlns:x0="http://www.w3.org/2001/XMLSchema">
  <time num_steps="2000" output_inc="2000" time_step="0.02">
    <schedule_function><piecewise_linear>
      <OrderedPair x="0" y="0"/><OrderedPair x="39" y="1"/><OrderedPair x="10000000" y="1"/>
    </piecewise_linear></schedule_function>
  </time>
  <nodes>
    <field field_name="electric_scalar_potential" integrator="static" solution_group="1">
      <dof_labels><String value="Psi"/></dof_labels>
      <kinematic_BC dof="1" node_ID="1" type="fixed"/>
      <kinematic_BC dof="1" node_ID="2" type="fixed"/>
    </field>
    <field field_name="displacement" integrator="static" solution_group="2">
      <dof_labels><String value="D_X"/><String value="D_Y"/></dof_labels>
      <kinematic_BC dof="2" node_ID="1" type="fixed"/>
      <kinematic_BC dof="1" node_ID="3" schedule="1" type="u" value="{g3:+.4f}"/>
      <kinematic_BC dof="1" node_ID="4" schedule="1" type="u" value="{g4:+.4f}"/>
    </field>
  </nodes>
  <element_list>
    <diffusion field_name="electric_scalar_potential"><quadrilateral/>
      <diffusion_element_block><block_ID_list><String value="1"/></block_ID_list>
        <diffusion_material><linear_dielectric_material epsilon="1.0"/></diffusion_material>
      </diffusion_element_block></diffusion>
    <updated_lagrangian_Q1P0_surface field_name="displacement" mass_type="lumped_mass" epsilon="1.0">
      <quadrilateral/><solid_element_nodal_output displacements="1" stress="1"/>
      <large_strain_element_block><block_ID_list><String value="1"/></block_ID_list>
        <large_strain_material_2D><RG_split_general constraint_2D="plane_strain" density="1.0">
          <rg_eq_potential><neo-hookean kappa="1000.67" mu="1.0"/></rg_eq_potential>
        </RG_split_general></large_strain_material_2D></large_strain_element_block>
      <surface_tension side_set_ID="1" gamma="0.0" t_0="20.0"/>
    </updated_lagrangian_Q1P0_surface>
  </element_list>
  <linear_solver><MUMPS_matrix/></linear_solver>
  <nonlinear_solver abs_tolerance="1.0e-7" check_LHS_perturbation="1.0E-8" check_code="no_check"
    divergence_tolerance="1.0e+10" max_iterations="100" rel_tolerance="1.0e-9"><MUMPS_matrix/></nonlinear_solver>
  <solver_phases max_loops="5">
    <solver_phase iterations="1" pass_iterations="1" solver="1"/>
    <solver_phase iterations="50" pass_iterations="1" solver="2"/>
  </solver_phases>
</tahoe>
"""


def stage2_xml(lam, gbar, nsteps=None):
    g3, g4 = (1.0 - lam) * GRIP, -(1.0 - lam) * GRIP
    gamma = 4.0 * gbar
    ns = nsteps if nsteps is not None else NSTEPS
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<tahoe author="sam" compute_IC="false" geometry_file="{MESH}"
    output_format="ExodusII" restart_file="stage1_l{lam_tag(lam)}.rs2000of2000"
    title="stage2 lam={lam} gbar={gbar}" xmlns:x0="http://www.w3.org/2001/XMLSchema">
  <time num_steps="{ns}" output_inc="50" time_step="{DT}">
    <schedule_function><piecewise_linear>
      <OrderedPair x="0" y="1"/><OrderedPair x="10000000" y="1"/>
    </piecewise_linear></schedule_function>
    <schedule_function><piecewise_linear>
      <OrderedPair x="0" y="0"/><OrderedPair x="{T_V0}" y="0"/>
      <OrderedPair x="{T_V1}" y="1"/><OrderedPair x="10000000" y="1"/>
    </piecewise_linear></schedule_function>
  </time>
  <nodes>
    <field field_name="electric_scalar_potential" integrator="static" solution_group="1">
      <dof_labels><String value="Psi"/></dof_labels>
      <kinematic_BC dof="1" node_ID="1" type="fixed"/>
      <kinematic_BC dof="1" node_ID="2" schedule="2" type="u" value="{PHI_MAX}"/>
    </field>
    <field field_name="displacement" integrator="central_difference" solution_group="2">
      <dof_labels><String value="D_X"/><String value="D_Y"/></dof_labels>
      <kinematic_BC dof="2" node_ID="1" type="fixed"/>
      <kinematic_BC dof="1" node_ID="3" schedule="1" type="u" value="{g3:+.4f}"/>
      <kinematic_BC dof="1" node_ID="4" schedule="1" type="u" value="{g4:+.4f}"/>
    </field>
  </nodes>
  <element_list>
    <diffusion field_name="electric_scalar_potential"><quadrilateral/>
      <diffusion_element_nodal_output displacement="1" coordinates="1"/>
      <diffusion_element_block><block_ID_list><String value="1"/></block_ID_list>
        <diffusion_material><linear_dielectric_material epsilon="1.0"/></diffusion_material>
      </diffusion_element_block></diffusion>
    <updated_lagrangian_Q1P0_surface field_name="displacement" mass_type="lumped_mass" epsilon="1.0">
      <quadrilateral/><solid_element_nodal_output displacements="1" stress="1"/>
      <large_strain_element_block><block_ID_list><String value="1"/></block_ID_list>
        <large_strain_material_2D><RG_split_general constraint_2D="plane_strain" density="1.0">
          <rg_eq_potential><neo-hookean kappa="1000.67" mu="1.0"/></rg_eq_potential>
        </RG_split_general></large_strain_material_2D></large_strain_element_block>
      <surface_tension side_set_ID="1" gamma="{gamma}" t_0="{T_GAMMA}"/>
    </updated_lagrangian_Q1P0_surface>
  </element_list>
  <linear_solver><MUMPS_matrix/></linear_solver>
  <linear_solver><diagonal_matrix/></linear_solver>
  <solver_phases max_loops="1">
    <solver_phase iterations="1" pass_iterations="0" solver="1"/>
    <solver_phase iterations="1" pass_iterations="0" solver="2"/>
  </solver_phases>
</tahoe>
"""


def run(deck):
    t0 = time.time()
    p = subprocess.run([TAHOE, "-f", deck], cwd=HERE,
                       stdout=open(os.path.join(HERE, deck + ".log"), "w"),
                       stderr=subprocess.STDOUT)
    return p.returncode, time.time() - t0


def all_lams():
    return sorted(set([1.0] + LAMS))


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "help"

    if mode == "stage1":
        for lam in all_lams():
            f = f"stage1_l{lam_tag(lam)}.xml"
            open(os.path.join(HERE, f), "w").write(stage1_xml(lam))
            rc, dt = run(f)
            print(f"stage1 lam={lam}: rc={rc} ({dt:.0f}s)", flush=True)

    elif mode == "stage2-decks":
        n = 0
        for gb in GBARS:                       # gamma sweep at lam=1
            f = f"stage2_l100_g{int(gb*10):03d}.xml"
            open(os.path.join(HERE, f), "w").write(stage2_xml(1.0, gb)); n += 1
        for lam in LAMS:                       # pre-strain sweep at LAM_GBAR
            f = f"stage2_l{lam_tag(lam)}_g{int(LAM_GBAR*10):03d}.xml"
            open(os.path.join(HERE, f), "w").write(stage2_xml(lam, LAM_GBAR)); n += 1
        print(f"{n} stage2 decks written (mesh={MESH})", flush=True)

    elif mode == "run":
        rc, dt = run(sys.argv[2])
        print(f"{sys.argv[2]}: rc={rc} ({dt:.0f}s)", flush=True)

    elif mode == "time":
        open(os.path.join(HERE, "_time.xml"), "w").write(stage2_xml(1.0, 1.0, nsteps=NSTEPS))
        rc, dt = run("_time.xml")
        print(f"timing: {NSTEPS} steps rc={rc} -> {dt:.1f}s "
              f"({1000*dt/NSTEPS:.2f} ms/step; full {int(T_V1/DT)+2000} steps ~ "
              f"{dt/NSTEPS*(int(T_V1/DT)+2000)/60:.1f} min)", flush=True)
    else:
        print(__doc__)


if __name__ == "__main__":
    main()
