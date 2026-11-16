import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
import numpy as np

from analyze import analyze, read
from sample import sample


def run(script, out, **variables):
    out.mkdir(parents=True, exist_ok=True)
    if (out / "done").exists():
        return
    command = [args.lmp, "-in", script, "-log", str(out / "log.lammps"),
               "-screen", "none", "-nocite", "-var", "out", str(out),
               "-var", "dt", str(settings["dt"])]
    for key, value in variables.items():
        command += ["-var", key, str(value)]
    start = time.monotonic()
    with (out / "error.out").open("w") as stream:
        subprocess.run(command, stdout=stream, stderr=stream, check=True)
    if script in ("in.ti", "in.widom"):
        sample(out, settings, variables["T"], variables.get("lam"))
    elapsed = time.monotonic() - start
    (out / "done").write_text(f"{elapsed:.1f}\n")
    print(f"{out}: {elapsed:.1f} s", flush=True)


def calculate(T):
    folder = root / str(T)
    run("in.relax", folder, data="h2o_equil.data", T=T,
        eq=100 if args.test else 50000)
    grid = [0.0, 0.5, 1.0] if args.test else np.linspace(0, 1, 21).tolist()
    for attempt in range(1 if args.test else 3):
        trial = folder / str(attempt)
        trial.mkdir(exist_ok=True)
        (trial / "lambda.json").write_text(json.dumps(grid) + "\n")
        data = folder / "equil.data"
        curves = []
        for direction, values in (("off", grid[::-1]), ("on", grid)):
            curve = []
            for lam in values:
                out = trial / direction / f"{lam:.6f}"
                run("in.ti", out, data=data, T=T, lam=lam, eq=settings["eq"],
                    steps=settings["steps"], stride=settings["stride"])
                data = out / "final.data"
                curve.append((lam, read(out / "deriv.out")[:, 1].mean()))
            curves.append(np.array(sorted(curve)))
        gap = abs(np.trapezoid(curves[0][:, 1] - curves[1][:, 1], grid))
        print(f"T={T}: off/on difference = {gap:.4f} kcal/mol", flush=True)
        if gap <= 1 or args.test:
            break
        if attempt < 2:
            y = np.mean([c[:, 1] for c in curves], axis=0)
            score = np.maximum(abs(y[:-1]), abs(y[1:])) + abs(np.diff(y))
            extra = [(grid[i] + grid[i + 1]) / 2 for i in np.argsort(score)[-6:]]
            grid = sorted(grid + extra)
    run("in.widom", folder / "widom", data=folder / "equil.data", T=T,
        eq=settings["eq"], steps=settings["steps"], stride=settings["stride"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--lmp", default=shutil.which("lmp_serial") or shutil.which("lmp"))
    args = parser.parse_args()
    if not args.lmp:
        parser.error("LAMMPS не найден: укажите --lmp")
    os.chdir(Path(__file__).resolve().parent)
    root = Path("check_optimized" if args.test else "results")
    settings = dict(eq=100 if args.test else 10000, steps=500 if args.test else 50000,
                    block=50 if args.test else 500, stride=5 if args.test else 50,
                    dt=2, insertions=10 if args.test else 1000, seed=340509,
                    test=args.test, lmp=args.lmp)
    if root.exists():
        if not args.resume or json.loads((root / "settings.json").read_text()) != settings:
            parser.error(f"{root} уже существует; --resume продолжит расчёт с теми же параметрами")
    root.mkdir(exist_ok=True)
    (root / "settings.json").write_text(json.dumps(settings, indent=2) + "\n")
    (root / "pid").write_text(str(os.getpid()) + "\n")
    (root / "status.txt").write_text("running\n")
    start = time.monotonic()
    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            list(pool.map(calculate, (300, 350)))
        analyze(root)
        (root / "status.txt").write_text("finished\n")
    except Exception:
        (root / "status.txt").write_text("failed; see run.log and window log.lammps\n")
        raise
    finally:
        (root / "time.txt").write_text(f"{time.monotonic() - start:.1f}\n")
