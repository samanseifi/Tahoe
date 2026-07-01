#!/usr/bin/env python3
"""
Clean electromechanical FEA sweep for the pre-strain trend.

Two-stage workflow per pre-stretch lambda:
  Stage 1 (static): grip the L=80, H=4 strip to in-plane stretch lambda, write restart.
  Stage 2 (explicit central-difference): restart, ramp surface tension then voltage;
           the wrinkle onset (top-surface undulation blow-up) gives the critical
           nominal field  Etilde_c*sqrt(eps/mu) = Phi/H,  H=4  (mu=eps=1).

Grip displacement for stretch lambda (strip length L=80, half each end):
  node 3 (left)  = +(1-lambda)*40      node 4 (right) = -(1-lambda)*40
gamma = 4*gbar  (since gbar = gamma/(mu H), H=4).  Voltage ramped to Phi_max so
Etilde reaches ~PHI_MAX/4 (set high enough to catch the stabilized stretch cases).
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TAHOE = "/home/samanseifi/codes/tahoe/build/bin/tahoe"
H = 4.0
PHI_MAX = 40.0          # Etilde_max = PHI_MAX/H = 10
DT = 0.02
T_GAMMA = 40.0          # gamma fully ramped by t=40
T_V0 = 40.0             # voltage starts
T_V1 = 1240.0           # voltage at full (Etilde=10) -> rate matches the validated decks
NSTEPS = int(T_V1 / DT) + 2000     # a little past full ramp


def lam_tag(lam):
    return f"{int(round(lam*100)):03d}"


def stage1_xml(lam):
    g3 = (1.0 - lam) * 40.0
    g4 = -(1.0 - lam) * 40.0
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<tahoe author="sam" compute_IC="false" geometry_file="../meshes/bar_2D.geom"
    output_format="ExodusII" restart_output_inc="2000" title="stage1 lam={lam}"
    xmlns:x0="http://www.w3.org/2001/XMLSchema">
  <time num_steps="2000" output_inc="2000" time_step="{DT}">
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


def stage2_xml(lam, gbar):
    g3 = (1.0 - lam) * 40.0
    g4 = -(1.0 - lam) * 40.0
    gamma = 4.0 * gbar
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<tahoe author="sam" compute_IC="false" geometry_file="../meshes/bar_2D.geom"
    output_format="ExodusII" restart_file="stage1_l{lam_tag(lam)}.rs2000of2000"
    title="stage2 lam={lam} gbar={gbar}" xmlns:x0="http://www.w3.org/2001/XMLSchema">
  <time num_steps="{NSTEPS}" output_inc="50" time_step="{DT}">
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
    p = subprocess.run([TAHOE, "-f", deck], cwd=HERE,
                       stdout=open(os.path.join(HERE, deck + ".log"), "w"),
                       stderr=subprocess.STDOUT)
    return p.returncode


def main():
    # matrix: pre-strain sweep at gbar=1, plus a couple of cross-checks at gbar=2
    LAMS = [0.85, 0.90, 0.95, 1.05, 1.10]
    GBARS = [1.0]
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"

    if mode in ("stage1", "all"):
        for lam in LAMS:
            f = f"stage1_l{lam_tag(lam)}.xml"
            open(os.path.join(HERE, f), "w").write(stage1_xml(lam))
            rc = run(f)
            print(f"stage1 lam={lam}: rc={rc}", flush=True)

    if mode in ("stage2", "all"):
        # write all stage2 decks; run them (caller may background)
        for lam in LAMS:
            for gb in GBARS:
                f = f"stage2_l{lam_tag(lam)}_g{int(gb*10):02d}.xml"
                open(os.path.join(HERE, f), "w").write(stage2_xml(lam, gb))
        print("stage2 decks written", flush=True)


if __name__ == "__main__":
    main()
