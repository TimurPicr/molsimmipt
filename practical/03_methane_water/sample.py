from itertools import islice, product

import numpy as np

R = 0.00198720425864083
EPS = np.sqrt(0.185207 * 0.294)
SIGMA = 3.44445
RC = 14.0


def frames(path):
    with open(path) as stream:
        while stream.readline():
            step = int(stream.readline())
            stream.readline()
            n = int(stream.readline())
            stream.readline()
            bounds = np.array([list(map(float, stream.readline().split())) for _ in range(3)])
            stream.readline()
            atoms = np.loadtxt(list(islice(stream, n)), ndmin=2)
            yield step, bounds[:, 0], bounds[:, 1] - bounds[:, 0], atoms


def energy(points, oxygen, box, lam=1.0, derivative=False):
    d = np.asarray(points).reshape(-1, 1, 3) - oxygen
    d -= box * np.rint(d / box)
    result = np.zeros(len(d))
    for shift in product((-1, 0, 1), repeat=3):
        r2 = np.sum((d + np.array(shift) * box) ** 2, axis=2)
        mask = r2 < RC ** 2
        a = np.maximum((r2 / SIGMA ** 2) ** 3 + 0.5 * (1 - lam) ** 2, 1e-30)
        inv = 1 / a
        u = inv ** 2 - inv
        if derivative:
            u = 2 * lam * u + lam ** 2 * (1 - lam) * (2 * inv ** 3 - inv ** 2)
        else:
            u *= lam ** 2
        result += 4 * EPS * np.where(mask, u, 0).sum(axis=1)
    tail = 16 * np.pi / 3 * len(oxygen) / np.prod(box) * EPS * SIGMA ** 3
    tail *= (SIGMA / RC) ** 9 / 3 - (SIGMA / RC) ** 3
    return result + tail * (2 * lam if derivative else lam ** 2)


def sample(out, settings, T, lam=None):
    rng = np.random.default_rng(settings["seed"] + T)
    values = []
    for step, lo, box, atoms in frames(out / "positions.dump"):
        if step == 0:
            continue
        oxygen = atoms[atoms[:, 1] == 1, 2:5]
        if lam is not None:
            methane = atoms[atoms[:, 1] == 3, 2:5]
            values.append(energy(methane, oxygen, box, lam, derivative=True)[0])
        else:
            total = 0.0
            for first in range(0, settings["insertions"], 128):
                n = min(128, settings["insertions"] - first)
                points = lo + rng.random((n, 3)) * box
                total += np.exp(-energy(points, oxygen, box) / (R * T)).sum()
            values.append(total / settings["insertions"])
    values = np.array(values)
    if lam is not None:
        count = settings["block"] // settings["stride"]
        means = values.reshape(-1, count).mean(axis=1)
        data = np.c_[np.arange(1, len(means) + 1) * settings["block"], means]
        np.savetxt(out / "deriv.out", data, header="step dU_dlambda")
    else:
        mu = np.zeros_like(values)
        positive = values > 0
        mu[positive] = -R * T * np.log(values[positive])
        data = np.c_[np.arange(1, len(values) + 1) * settings["stride"], mu, values]
        np.savetxt(out / "widom.out", data, header="step mu mean_weight")
