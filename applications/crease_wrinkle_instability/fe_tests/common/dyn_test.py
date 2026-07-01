"""Dynamic-effects test: lam=1, gbar=1, coarse mesh, 3 ramp rates.
If onset field is rate-independent -> quasi-static (dynamics NOT the cause)."""
import os, subprocess, time, numpy as np
from netCDF4 import Dataset
HERE=os.path.dirname(os.path.abspath(__file__)); TAHOE="/home/samanseifi/codes/tahoe/build/bin/tahoe"
MESH="../../meshes/bar_2D.geom"; DT=0.02; PHI=40.0; H=4.0
def stage1():
    return f'''<?xml version="1.0"?><tahoe author="s" compute_IC="false" geometry_file="{MESH}" output_format="ExodusII" restart_output_inc="2000" title="s1dyn" xmlns:x0="http://www.w3.org/2001/XMLSchema">
<time num_steps="2000" output_inc="2000" time_step="0.02"><schedule_function><piecewise_linear><OrderedPair x="0" y="0"/><OrderedPair x="39" y="1"/><OrderedPair x="1e7" y="1"/></piecewise_linear></schedule_function></time>
<nodes><field field_name="electric_scalar_potential" integrator="static" solution_group="1"><dof_labels><String value="Psi"/></dof_labels><kinematic_BC dof="1" node_ID="1" type="fixed"/><kinematic_BC dof="1" node_ID="2" type="fixed"/></field>
<field field_name="displacement" integrator="static" solution_group="2"><dof_labels><String value="D_X"/><String value="D_Y"/></dof_labels><kinematic_BC dof="2" node_ID="1" type="fixed"/><kinematic_BC dof="1" node_ID="3" type="fixed"/><kinematic_BC dof="1" node_ID="4" type="fixed"/></field></nodes>
<element_list><diffusion field_name="electric_scalar_potential"><quadrilateral/><diffusion_element_block><block_ID_list><String value="1"/></block_ID_list><diffusion_material><linear_dielectric_material epsilon="1.0"/></diffusion_material></diffusion_element_block></diffusion>
<updated_lagrangian_Q1P0_surface field_name="displacement" mass_type="lumped_mass" epsilon="1.0"><quadrilateral/><solid_element_nodal_output displacements="1" stress="1"/><large_strain_element_block><block_ID_list><String value="1"/></block_ID_list><large_strain_material_2D><RG_split_general constraint_2D="plane_strain" density="1.0"><rg_eq_potential><neo-hookean kappa="1000.67" mu="1.0"/></rg_eq_potential></RG_split_general></large_strain_material_2D></large_strain_element_block><surface_tension side_set_ID="1" gamma="0.0" t_0="20.0"/></updated_lagrangian_Q1P0_surface></element_list>
<linear_solver><MUMPS_matrix/></linear_solver><nonlinear_solver abs_tolerance="1e-7" check_LHS_perturbation="1e-8" check_code="no_check" divergence_tolerance="1e10" max_iterations="100" rel_tolerance="1e-9"><MUMPS_matrix/></nonlinear_solver>
<solver_phases max_loops="5"><solver_phase iterations="1" pass_iterations="1" solver="1"/><solver_phase iterations="50" pass_iterations="1" solver="2"/></solver_phases></tahoe>'''
def stage2(tv1, nsteps, rsname):
    return f'''<?xml version="1.0"?><tahoe author="s" compute_IC="false" geometry_file="{MESH}" output_format="ExodusII" restart_file="{rsname}" title="s2dyn" xmlns:x0="http://www.w3.org/2001/XMLSchema">
<time num_steps="{nsteps}" output_inc="50" time_step="{DT}"><schedule_function><piecewise_linear><OrderedPair x="0" y="1"/><OrderedPair x="1e7" y="1"/></piecewise_linear></schedule_function><schedule_function><piecewise_linear><OrderedPair x="0" y="0"/><OrderedPair x="40" y="0"/><OrderedPair x="{tv1}" y="1"/><OrderedPair x="1e7" y="1"/></piecewise_linear></schedule_function></time>
<nodes><field field_name="electric_scalar_potential" integrator="static" solution_group="1"><dof_labels><String value="Psi"/></dof_labels><kinematic_BC dof="1" node_ID="1" type="fixed"/><kinematic_BC dof="1" node_ID="2" schedule="2" type="u" value="{PHI}"/></field>
<field field_name="displacement" integrator="central_difference" solution_group="2"><dof_labels><String value="D_X"/><String value="D_Y"/></dof_labels><kinematic_BC dof="2" node_ID="1" type="fixed"/><kinematic_BC dof="1" node_ID="3" type="fixed"/><kinematic_BC dof="1" node_ID="4" type="fixed"/></field></nodes>
<element_list><diffusion field_name="electric_scalar_potential"><quadrilateral/><diffusion_element_nodal_output displacement="1" coordinates="1"/><diffusion_element_block><block_ID_list><String value="1"/></block_ID_list><diffusion_material><linear_dielectric_material epsilon="1.0"/></diffusion_material></diffusion_element_block></diffusion>
<updated_lagrangian_Q1P0_surface field_name="displacement" mass_type="lumped_mass" epsilon="1.0"><quadrilateral/><solid_element_nodal_output displacements="1" stress="1"/><large_strain_element_block><block_ID_list><String value="1"/></block_ID_list><large_strain_material_2D><RG_split_general constraint_2D="plane_strain" density="1.0"><rg_eq_potential><neo-hookean kappa="1000.67" mu="1.0"/></rg_eq_potential></RG_split_general></large_strain_material_2D></large_strain_element_block><surface_tension side_set_ID="1" gamma="4.0" t_0="40.0"/></updated_lagrangian_Q1P0_surface></element_list>
<linear_solver><MUMPS_matrix/></linear_solver><linear_solver><diagonal_matrix/></linear_solver>
<solver_phases max_loops="1"><solver_phase iterations="1" pass_iterations="0" solver="1"/><solver_phase iterations="1" pass_iterations="0" solver="2"/></solver_phases></tahoe>'''
def run(deck):
    subprocess.run([TAHOE,"-f",deck],cwd=HERE,stdout=open(deck+".log","w"),stderr=subprocess.STDOUT)
open("s1dyn.xml","w").write(stage1()); run("s1dyn.xml")
import glob; rs=[x for x in glob.glob("s1dyn.rs*") if "." not in x.split("rs")[-1]][0]
print("restart:",rs)
def onset(f,tv1):
    ds=Dataset(f); cx=np.array(ds["coordx"][:]);cy=np.array(ds["coordy"][:]);t=np.array(ds["time_whole"][:])
    DX=np.array(ds["vals_nod_var1"][:]);DY=np.array(ds["vals_nod_var2"][:])
    top=np.where(np.isclose(cy,cy.max(),atol=1e-6))[0];x0=cx[top];Lx=x0.max()-x0.min()
    I=(x0>x0.min()+0.12*Lx)&(x0<x0.max()-0.12*Lx)
    amp=np.array([np.ptp((cy[top][I]+DY[k,top][I])-np.polyval(np.polyfit(x0[I]+DX[k,top][I],cy[top][I]+DY[k,top][I],2),x0[I]+DX[k,top][I])) for k in range(len(t))])
    E=np.clip((t-40.0)/(tv1-40.0),0,1)*10.0; on=np.where((amp>0.3)&(E>0.05))[0]
    return E[on[0]] if len(on) else np.nan
for tag,tv1 in (("fast",620.0),("base",1240.0),("slow",2480.0)):
    ns=int(tv1/DT)+2000; d=f"s2dyn_{tag}.xml"; open(d,"w").write(stage2(tv1,ns,rs))
    t0=time.time(); run(d); Eon=onset(d.replace(".xml",".io1.exo"),tv1)
    print(f"  {tag:4s} ramp T_V1={tv1:5.0f} (rate {10/(tv1-40):.4f}/t): Eon={Eon:.3f}  ({time.time()-t0:.0f}s)")
