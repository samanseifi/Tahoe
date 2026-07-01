#!/usr/bin/env bash
# Extended refined-mesh EM sweep driver.
#  stage1 (static, DT=0.02) for each needed lambda  ->  restart
#  stage2 (explicit, DT=0.005, capped at Etilde~7) launched in parallel
set -e
cd "$(dirname "$0")"
export MESH=../meshes/bar_2D_fineY.geom     # 80x16 (4x through-thickness refinement)
export DT=0.005
export STEPS=176000                          # Etilde ramps to 7 (covers all FE onsets)
PY=python3

echo "=== stage1 (static restarts) ==="
# write+run stage1 decks via the module (tags: 100,080,090,110,120)
for lam in 1.0 0.8 0.9 1.1 1.2; do
  tag=$($PY -c "import run_sweep2 as r;print(r.lam_tag($lam))")
  $PY -c "import run_sweep2 as r;open('stage1_l${tag}.xml','w').write(r.stage1_xml($lam))"
  $PY run_sweep2.py run stage1_l${tag}.xml
done

echo "=== write stage2 decks ==="
# gamma sweep at lam=1
for g in 0.5 1.0 2.0 5.0 10.0 20.0; do
  gg=$($PY -c "print(f'{int($g*10):03d}')")
  $PY -c "import run_sweep2 as r;open('stage2_l100_g${gg}.xml','w').write(r.stage2_xml(1.0,$g))"
done
# lambda sweep at gbar=2
for lam in 0.8 0.9 1.1 1.2; do
  tag=$($PY -c "import run_sweep2 as r;print(r.lam_tag($lam))")
  $PY -c "import run_sweep2 as r;open('stage2_l${tag}_g020.xml','w').write(r.stage2_xml($lam,2.0))"
done

echo "=== launch stage2 in parallel ==="
ls stage2_l100_g005.xml stage2_l100_g010.xml stage2_l100_g020.xml \
   stage2_l100_g050.xml stage2_l100_g100.xml stage2_l100_g200.xml \
   stage2_l080_g020.xml stage2_l090_g020.xml stage2_l110_g020.xml stage2_l120_g020.xml \
   | xargs -P 10 -I{} $PY run_sweep2.py run {}
echo "=== ALL DONE ==="
