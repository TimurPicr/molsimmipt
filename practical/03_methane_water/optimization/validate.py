from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sample import energy, frames

root = Path('optimization')


def run(name, text, dt=2):
    path = root / f'{name}.in'
    path.write_text(text)
    subprocess.run(['lmp_serial', '-in', str(path), '-log', str(root / f'{name}.log'),
                    '-screen', 'none', '-nocite', '-var', 'data', 'h2o_equil.data',
                    '-var', 'dt', str(dt)], check=True)


checks = []
for style, common, coeff in (
    ('old', 'optimization/in.common.original', ''),
    ('fast', 'in.common', 'lj/cut/soft '),
):
    for lam in (0, 0.2, 0.5, 1):
        name = f'{style}_{lam}'
        run(name, f'''include {common}
unfix myshk
pair_coeff 1 3 {coeff}0.23334707626194934 3.44445 {lam}
run 0
write_dump all custom optimization/{name}.dump id type x y z fx fy fz modify sort id format float %.15g
print "$(pe:%.16g)" file optimization/{name}.energy screen no
''')

step, lo, box, atoms = next(frames(root / 'fast_1.dump'))
oxygen = atoms[atoms[:, 1] == 1, 2:5]
methane = atoms[atoms[:, 1] == 3, 2:5]
base = float((root / 'fast_0.energy').read_text())
for lam in (0, 0.2, 0.5, 1):
    actual = float((root / f'fast_{lam}.energy').read_text()) - base
    predicted = energy(methane, oxygen, box, lam)[0]
    assert abs(actual - predicted) < 1e-7, (lam, actual, predicted)
    h = 1e-5
    a, b = max(0, lam - h), min(1, lam + h)
    finite = (energy(methane, oxygen, box, b) - energy(methane, oxygen, box, a))[0] / (b - a)
    exact = energy(methane, oxygen, box, lam, derivative=True)[0]
    assert abs(exact - finite) < 0.002, (lam, exact, finite)
    old = next(frames(root / f'old_{lam}.dump'))[3]
    fast = next(frames(root / f'fast_{lam}.dump'))[3]
    relative = np.linalg.norm(old[:, 5:8] - fast[:, 5:8]) / np.linalg.norm(old[:, 5:8])
    assert relative < 0.002, relative
    checks.append(f'lambda={lam}: energy error={actual-predicted:.3g}; derivative error={exact-finite:.3g}; force relative RMS={relative:.3g}')

rng = np.random.default_rng(340509)
for i, point in enumerate(lo + rng.random((5, 3)) * box):
    x, y, z = point
    run(f'insert_{i}', f'''include in.common
unfix myshk
set group CH4 x {x} y {y} z {z}
pair_coeff 1 3 lj/cut/soft 0.23334707626194934 3.44445 0.0
run 0
variable zero equal $(pe:%.16g)
pair_coeff 1 3 lj/cut/soft 0.23334707626194934 3.44445 1.0
run 0
print "$(pe-v_zero:%.16g)" file optimization/insert_{i}.energy screen no
''')
    actual = float((root / f'insert_{i}.energy').read_text())
    predicted = energy(point, oxygen, box)[0]
    assert np.isclose(actual, predicted, rtol=1e-10, atol=1e-7), (actual, predicted)
    checks.append(f'insertion {i}: energy error={actual-predicted:.3g} kcal/mol')

for dt in (1, 2):
    run(f'nve_{dt}', f'''include in.common
fix integ all nvt temp 300 300 100
run 100
unfix integ
fix integ all nve
reset_timestep 0
variable total equal etotal
fix avg all ave/time 10 1 10 v_total file optimization/nve_{dt}.out
run {2000 // dt}
''', dt)
    data = np.loadtxt(root / f'nve_{dt}.out')
    drift = np.polyfit(data[:, 0] * dt / 1000, data[:, 1], 1)[0]
    assert abs(drift) < 0.5, (dt, drift)
    checks.append(f'dt={dt} fs: NVE energy drift={drift:.4g} kcal/mol/ps')
(root / 'checks.txt').write_text('\n'.join(checks) + '\n')
print('\n'.join(checks))
