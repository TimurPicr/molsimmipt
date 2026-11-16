import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

R = 0.00198720425864083


def read(path):
    data = np.loadtxt(path, comments="#", ndmin=2)
    data = data[data[:, 0] > 0]
    if not len(data) or not np.isfinite(data).all():
        raise ValueError(f"Нет корректных данных: {path}")
    return data


def blocks(values, n=10):
    return np.array([a.mean() for a in np.array_split(values, min(n, len(values)))])


def error(values):
    return max(np.std(b, ddof=1) / np.sqrt(len(b))
               for b in (blocks(values, 5), blocks(values, 10)))


def analyze(root):
    root = Path(root)
    settings = json.loads((root / "settings.json").read_text())
    rows, energies, widom = [], [], []
    text = ["Проверочный прогон; физические выводы делать нельзя.\n" if settings["test"]
            else "Энергия сольватации метана, ккал/моль; ошибки — 1 SEM.\n"]
    for T in (300, 350):
        folder = root / str(T)
        trial = max((p for p in folder.iterdir() if p.name.isdigit()), key=lambda p: int(p.name))
        grid = np.array(json.loads((trial / "lambda.json").read_text()))
        weights = np.r_[np.diff(grid)[0], grid[2:] - grid[:-2], np.diff(grid)[-1]] / 2
        estimates = []
        fig, ax = plt.subplots()
        for direction in ("off", "on"):
            samples = [read(trial / direction / f"{lam:.6f}" / "deriv.out")[:, 1] for lam in grid]
            means = np.array([a.mean() for a in samples])
            errors = np.array([error(a) for a in samples])
            g = weights @ means
            se = np.linalg.norm(weights * errors)
            estimates.append((g, se))
            rows.append((T, direction, g, se))
            ax.errorbar(grid, means, yerr=errors, marker=".", label=direction)
            np.savetxt(trial / f"{direction}.csv", np.c_[grid, means, errors],
                       delimiter=",", header="lambda,dU_dlambda,sem", comments="")
        ax.set(xlabel="lambda", ylabel="<dU/dlambda>, kcal/mol", title=f"{T} K")
        ax.legend()
        fig.tight_layout()
        fig.savefig(folder / "ti.png", dpi=160)
        plt.close(fig)
        gap = abs(estimates[0][0] - estimates[1][0])
        g = np.mean([a[0] for a in estimates])
        se = np.linalg.norm([a[1] for a in estimates]) / 2
        energies.append((g, se))
        text.append(f"{T} K: off {estimates[0][0]:.3f} ± {estimates[0][1]:.3f}; "
                    f"on {estimates[1][0]:.3f} ± {estimates[1][1]:.3f}; разность {gap:.3f}.")
        if gap > 1:
            text.append("TI не сошлось: требуется дальнейшее уточнение сетки и увеличение выборки.")
        data = read(folder / "widom/widom.out")
        values = data[:, 2]
        if np.any(values < 0) or values.mean() <= 0:
            raise ValueError("Не накоплены веса Уидома")
        positive = values > 0
        if not np.allclose(data[positive, 1], -R * T * np.log(values[positive]), rtol=2e-5, atol=2e-5):
            raise ValueError("Проверьте формат f_w[2] в установленном LAMMPS")
        mu = -R * T * np.log(values.mean())
        widom.append(mu)
        se_w = R * T * error(values / values.mean())
        rows.append((T, "widom", mu, se_w))
        text.append(f"Уидом: {mu:.3f} ± {se_w:.3f}; отличие от среднего TI {mu - g:+.3f}.")
        fig, ax = plt.subplots()
        for direction in ("off", "on"):
            samples = [read(trial / direction / f"{lam:.6f}" / "deriv.out")[:, 1] for lam in grid]
            cumulative = np.array([np.cumsum(a) / np.arange(1, len(a) + 1) for a in samples])
            ax.plot(np.arange(1, cumulative.shape[1] + 1) * settings["block"] * settings.get("dt", 1) / 1000,
                    weights @ cumulative, label=f"TI {direction}")
        cumulative = np.cumsum(values) / np.arange(1, len(values) + 1)
        good = cumulative > 0
        ax.plot(np.arange(1, len(values) + 1)[good] * settings.get("stride", 100) * settings.get("dt", 1) / 1000,
                -R * T * np.log(cumulative[good]), label="Widom")
        ax.set(xlabel="sampling time, ps", ylabel="mu_ex, kcal/mol", title=f"{T} K")
        ax.legend()
        fig.tight_layout()
        fig.savefig(folder / "convergence.png", dpi=160)
        plt.close(fig)
    (g1, e1), (g2, e2) = energies
    entropy = -(g2 - g1) / 50
    entropy_error = np.hypot(e1, e2) / 50
    enthalpy = (g1 + g2) / 2 + 325 * entropy
    enthalpy_error = np.hypot(7 * e1, 6 * e2)
    text += [f"325 K: ΔS = {entropy * 1000:.2f} ± {entropy_error * 1000:.2f} кал/(моль·К); "
             f"ΔH = {enthalpy:.3f} ± {enthalpy_error:.3f} ккал/моль.",
             "Погрешности статистические, без ошибки сетки и модели; независимость окон приближённая."]
    kh = 0.0014 * np.exp(1600 * (1 / 300 - 1 / 298.15))
    experiment = -R * 300 * np.log(8.314462618 * 300 * kh * 996.556 / 100000)
    text += [f"Эксперимент при 300 K: ≈{experiment:.2f} ккал/моль по константе Генри NIST; "
             f"TI − эксперимент = {g1 - experiment:+.3f}; "
             f"Уидом − эксперимент = {widom[0] - experiment:+.3f} ккал/моль.",
             "Для 350 K экспериментальное значение не экстраполировано за пределы применимости данных."]
    with (root / "summary.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(("T", "method", "mu_ex_kcal_mol", "sem_kcal_mol"))
        writer.writerows(rows)
    (root / "summary.txt").write_text("\n".join(text) + "\n")
    if not settings["test"]:
        report = Path(__file__).with_name("REPORT.md")
        intro = report.read_text().split("\nРезультаты расчёта\n")[0]
        report.write_text(intro + "\nРезультаты расчёта\n\n" + "\n".join(text) + "\n")
    print("\n".join(text))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", nargs="?", default="results")
    analyze(parser.parse_args().folder)
