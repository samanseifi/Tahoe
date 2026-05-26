#!/usr/bin/env bash
# sweep_prestretch_phase_diagram.sh
#
# Two-stage pre-stretch + voltage sweep over (eps_pre, gamma) to map
# the crease/wrinkle V_crit phase diagram for paper #1.
#
# Workflow per (eps_pre, gamma) point:
#   Stage 1: static mech pre-stretch with gamma=0 (no surface tension),
#            both fields declared so Stage 2 can restart cleanly.
#            dt=0.01 for Newton robustness across eps_pre range.
#
#   Stage 2: staggered explicit (central-difference mech + diffusion Psi),
#            restart from Stage 1, ramp gamma via t_0 then ramp V.
#            dt=0.02 (CFL-safe), V_max chosen ~1.5x analytic V_crease.
#
# Per-point wall: ~6-10 min.  Full sweep over ~20 points: ~3 hr.
#
# Usage:
#   ./sweep_prestretch_phase_diagram.sh                  # full sweep
#   ./sweep_prestretch_phase_diagram.sh 0.10 2.0         # one point
#
# Outputs (per point): swp_eps{NN}_g{GG}.{stage1,stage2}.{out,log,io*.exo,rs*}
# Plus swp_summary.csv with columns: eps_pre, gamma, V_target, last_t, last_V,
# d_y_spread_max, status (completed / element_inversion / max_iter)

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
TAHOE=${TAHOE:-$REPO_ROOT/build/bin/tahoe}
# Generated XMLs are run from $RUNROOT/<tag> (CWD), so they need a path that
# resolves to the meshes/ dir regardless of CWD depth — use absolute.
MESH="$(realpath "$SCRIPT_DIR/../meshes/bar_2D.geom")"
H=4.0
L=80.0
RUNROOT="$SCRIPT_DIR/sweep_runs"
SUMMARY=$RUNROOT/swp_summary.csv
mkdir -p $RUNROOT

# Sweep grid - edit these for different sweeps.
EPS_LIST=(0.00 0.05 0.10 0.15 0.20 0.25 0.30)   # 7 values (added 25, 30)
GAM_LIST=(0.0 1.0 2.0 4.0)                       # 4 values (gamma_bar = 0, 0.25, 0.5, 1)

# Override with command-line: one point
if [ $# -eq 2 ]; then
    EPS_LIST=("$1")
    GAM_LIST=("$2")
fi

# Init summary
[ -f $SUMMARY ] || echo "eps_pre,gamma,V_target,last_t,last_V,d_y_spread_max,status" > $SUMMARY

write_stage1_xml() {
    local eps=$1; local pre=$(awk -v e=$eps -v L=$L 'BEGIN{printf "%.4f", e*L/2}')
    local tag=$2
    cat > $RUNROOT/$tag/swp_${tag}.stage1.xml <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<tahoe author="sweep" compute_IC="false" geometry_file="$MESH" output_format="ExodusII"
    restart_output_inc="2000"
    title="Stage1 eps_pre=$eps (Q1P0_surface, gamma=0)" xmlns:x0="http://www.w3.org/2001/XMLSchema">
    <time num_steps="2000" output_inc="200" time_step="0.02">
        <schedule_function>
            <piecewise_linear>
                <OrderedPair x="0" y="0"/><OrderedPair x="20" y="1"/><OrderedPair x="10000000" y="1"/>
            </piecewise_linear>
        </schedule_function>
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
            <kinematic_BC dof="1" node_ID="3" schedule="1" type="u" value="-${pre}"/>
            <kinematic_BC dof="1" node_ID="4" schedule="1" type="u" value="+${pre}"/>
        </field>
    </nodes>
    <element_list>
        <diffusion field_name="electric_scalar_potential">
            <quadrilateral/>
            <diffusion_element_block>
                <block_ID_list><String value="1"/></block_ID_list>
                <diffusion_material><linear_dielectric_material epsilon="1.0"/></diffusion_material>
            </diffusion_element_block>
        </diffusion>
        <updated_lagrangian_Q1P0_surface field_name="displacement" mass_type="lumped_mass" epsilon="1.0">
            <quadrilateral/>
            <solid_element_nodal_output displacements="1" stress="1"/>
            <large_strain_element_block>
                <block_ID_list><String value="1"/></block_ID_list>
                <large_strain_material_2D>
                    <RG_split_general constraint_2D="plane_strain" density="1.0">
                        <rg_eq_potential><neo-hookean kappa="1000.67" mu="1.0"/></rg_eq_potential>
                    </RG_split_general>
                </large_strain_material_2D>
            </large_strain_element_block>
            <surface_tension side_set_ID="1" gamma="0.0" t_0="20.0"/>
        </updated_lagrangian_Q1P0_surface>
    </element_list>
    <linear_solver><MUMPS_matrix/></linear_solver>
    <nonlinear_solver abs_tolerance="1.0e-7" check_LHS_perturbation="1.0E-8" check_code="no_check"
        divergence_tolerance="1.0e+10" max_iterations="100" rel_tolerance="1.0e-9">
        <MUMPS_matrix/>
    </nonlinear_solver>
    <solver_phases max_loops="5">
        <solver_phase iterations="1"  pass_iterations="1" solver="1"/>
        <solver_phase iterations="50" pass_iterations="1" solver="2"/>
    </solver_phases>
</tahoe>
EOF
}

write_stage2_xml() {
    local eps=$1; local gam=$2; local pre=$(awk -v e=$eps -v L=$L 'BEGIN{printf "%.4f", e*L/2}')
    local tag=$3
    # V_target heuristic: analytic V_crease + headroom.  V_c ~ H_eff * (1.03 + 1.88*sqrt(gamma/(H_eff))).
    local v_target=$(awk -v e=$eps -v g=$gam -v H=$H 'BEGIN{
        heff = H/(1+e); gbar = g/heff; vc = heff*(1.03 + 1.88*sqrt(gbar));
        v = vc * 1.5; if (v < 4) v = 4; printf "%.2f", v;
    }')
    local t_end=$(awk -v v=$v_target 'BEGIN{printf "%.0f", 80 + v*50}')   # V_dot = 0.02
    local nsteps=$(awk -v t=$t_end 'BEGIN{printf "%d", (t-40)/0.02}')
    cat > $RUNROOT/$tag/swp_${tag}.stage2.xml <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<tahoe author="sweep" compute_IC="false" geometry_file="$MESH" output_format="ExodusII"
    restart_file="swp_${tag}.stage1.rs2000of2000"
    title="Stage2 eps=$eps gamma=$gam (explicit staggered)" xmlns:x0="http://www.w3.org/2001/XMLSchema">
    <time num_steps="$nsteps" output_inc="100" time_step="0.02">
        <schedule_function>
            <piecewise_linear><OrderedPair x="0" y="1"/><OrderedPair x="10000000" y="1"/></piecewise_linear>
        </schedule_function>
        <schedule_function>
            <piecewise_linear>
                <OrderedPair x="0" y="0"/><OrderedPair x="80" y="0"/>
                <OrderedPair x="$t_end" y="1"/><OrderedPair x="10000000" y="1"/>
            </piecewise_linear>
        </schedule_function>
    </time>
    <nodes>
        <field field_name="electric_scalar_potential" integrator="static" solution_group="1">
            <dof_labels><String value="Psi"/></dof_labels>
            <kinematic_BC dof="1" node_ID="1" type="fixed"/>
            <kinematic_BC dof="1" node_ID="2" schedule="2" type="u" value="${v_target}"/>
        </field>
        <field field_name="displacement" integrator="central_difference" solution_group="2">
            <dof_labels><String value="D_X"/><String value="D_Y"/></dof_labels>
            <kinematic_BC dof="2" node_ID="1" type="fixed"/>
            <kinematic_BC dof="1" node_ID="3" schedule="1" type="u" value="-${pre}"/>
            <kinematic_BC dof="1" node_ID="4" schedule="1" type="u" value="+${pre}"/>
        </field>
    </nodes>
    <element_list>
        <diffusion field_name="electric_scalar_potential">
            <quadrilateral/>
            <diffusion_element_block>
                <block_ID_list><String value="1"/></block_ID_list>
                <diffusion_material><linear_dielectric_material epsilon="1.0"/></diffusion_material>
            </diffusion_element_block>
        </diffusion>
        <updated_lagrangian_Q1P0_surface field_name="displacement" mass_type="lumped_mass" epsilon="1.0">
            <quadrilateral/>
            <solid_element_nodal_output displacements="1" stress="1"/>
            <large_strain_element_block>
                <block_ID_list><String value="1"/></block_ID_list>
                <large_strain_material_2D>
                    <RG_split_general constraint_2D="plane_strain" density="1.0">
                        <rg_eq_potential><neo-hookean kappa="1000.67" mu="1.0"/></rg_eq_potential>
                    </RG_split_general>
                </large_strain_material_2D>
            </large_strain_element_block>
            <surface_tension side_set_ID="1" gamma="${gam}" t_0="80.0"/>
        </updated_lagrangian_Q1P0_surface>
    </element_list>
    <linear_solver><MUMPS_matrix/></linear_solver>
    <linear_solver><diagonal_matrix/></linear_solver>
    <solver_phases max_loops="1">
        <solver_phase iterations="1" pass_iterations="0" solver="1"/>
        <solver_phase iterations="1" pass_iterations="0" solver="2"/>
    </solver_phases>
</tahoe>
EOF
    echo $v_target
}

analyze_stage2() {
    local tag=$1; local eps=$2; local gam=$3; local v_target=$4
    python3 - <<PY 2>/dev/null
import netCDF4 as nc, numpy as np, sys
try:
    ds = nc.Dataset("$RUNROOT/$tag/swp_${tag}.stage2.io1.exo")
    names = [b"".join(c).rstrip(b"\x00").decode() for c in ds.variables["name_nod_var"][:]]
    t = ds.variables["time_whole"][:]
    Y = ds.variables["coordy"][:]
    iy = names.index("D_Y")
    top = np.where(np.abs(Y - $H) < 1e-6)[0]
    last = ds.variables[f"vals_nod_var{iy+1}"][-1, top]
    spread_max = 0.0
    for k in range(len(t)):
        dy = ds.variables[f"vals_nod_var{iy+1}"][k, top]
        s = dy.max() - dy.min()
        if s > spread_max: spread_max = s
    last_t = float(t[-1])
    last_V = max(0, $v_target * (last_t-80)/($v_target*50))
    print(f"$eps,$gam,$v_target,{last_t:.2f},{last_V:.3f},{spread_max:.4e}")
except Exception as e:
    print(f"$eps,$gam,$v_target,nan,nan,nan", file=sys.stderr)
PY
}

# Main sweep loop
for eps in "${EPS_LIST[@]}"; do
    for gam in "${GAM_LIST[@]}"; do
        eps_pct=$(awk -v e=$eps 'BEGIN{printf "%02.0f", e*100}')
        gam_str=$(awk -v g=$gam 'BEGIN{printf "%02.0f", g*10}')
        tag="eps${eps_pct}_g${gam_str}"

        echo "=== ($eps, $gam) tag=$tag ==="
        mkdir -p $RUNROOT/$tag

        # Stage 1
        if [ ! -f $RUNROOT/$tag/swp_${tag}.stage1.rs2000of2000 ]; then
            write_stage1_xml $eps $tag
            echo "  Stage 1..."
            (cd $RUNROOT/$tag && $TAHOE -f swp_${tag}.stage1.xml > swp_${tag}.stage1.log 2>&1) || {
                echo "  Stage 1 FAILED for ($eps, $gam)"
                echo "${eps},${gam},nan,nan,nan,nan,stage1_failed" >> $SUMMARY
                continue
            }
            # Confirm restart file exists (tahoe can exit 0 even after element inversion)
            if [ ! -f $RUNROOT/$tag/swp_${tag}.stage1.rs2000of2000 ]; then
                echo "  Stage 1 INVERTED for ($eps, $gam) — no restart written"
                echo "${eps},${gam},nan,nan,nan,nan,stage1_inverted" >> $SUMMARY
                continue
            fi
        else
            echo "  Stage 1 cached"
        fi

        # Stage 2
        echo "  Stage 2..."
        v_target=$(write_stage2_xml $eps $gam $tag)
        if (cd $RUNROOT/$tag && $TAHOE -f swp_${tag}.stage2.xml > swp_${tag}.stage2.log 2>&1); then
            status="completed"
        else
            if grep -q "zero or negative jacobian" $RUNROOT/$tag/swp_${tag}.stage2.log 2>/dev/null; then
                status="element_inversion"
            else
                status="other_fail"
            fi
        fi

        # Extract last point + max spread
        line=$(analyze_stage2 $tag $eps $gam $v_target)
        echo "${line},${status}" >> $SUMMARY
        echo "  done: ${line},${status}"
    done
done

echo
echo "Sweep done.  Summary: $SUMMARY"
column -s, -t $SUMMARY | head -50
