"""Python port of the HIST-PPM-V3 space-time ETAS log-likelihood.

The conditional intensity is the one coded in func17 of hist-etas5pa.f,
hist-etas-mk.f and st-etas.f:

    lambda(t,x,y) = mu(x,y)
        + sum_{t_j<t} K_j / (t-t_j+c)^p_j
                      * ( Q_j(dx e^{-a_j M_j}, dy e^{-a_j M_j}) + d )^(-q_j)

with Q_j the unit-determinant anisotropy form of event j, M_j = mag - Mc, and
a_j the *halved* alpha that the programs use internally (a4 = aa4/2).  x is
longitude scaled by cos((YMIN+TY/2) deg), exactly as the Fortran does, so
every spatial quantity here is in (scaled) degrees.

Two evaluation modes:

``fortran``  reproduces a given program's number to rounding, quirks included
             (first-event term, 8 or 128 boundary segments, the sign of the
             x or y offset on the left and bottom sides of the box, and for
             st-etas the absolute coordinates and the uncompensated
             precursory events).  Used only to validate this port.
``exact``    the intensity as defined, with the compensator taken over the
             target period [tstart, tz] and the region box, the boundary
             integral on a fine segmentation with the correct offsets.
             Used for the residual analysis.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

FORTRAN_PI = 3.14159  # the constant the programs use for the cos-latitude factor


# ---------------------------------------------------------------- data
@dataclass
class Catalog:
    x: np.ndarray       # scaled longitude offset from XMIN (deg * cosfac)
    y: np.ndarray       # latitude offset from YMIN (deg)
    mag: np.ndarray     # magnitude minus Mc
    t: np.ndarray       # days from catalogue origin
    qa: np.ndarray      # anisotropy form  Q = qa x^2 - 2 qb x y + qc y^2
    qb: np.ndarray
    qc: np.ndarray
    tx: float           # scaled box width
    ty: float
    tz: float
    tstart: float
    xmg0: float
    xmin: float
    ymin: float
    cosfac: float
    # Delaunay nodes (events first, then the extra hull nodes) and triangles
    node_x: np.ndarray | None = None
    node_y: np.ndarray | None = None
    tri: np.ndarray | None = None       # (ndelo, 3) zero-based node indices
    tri_area: np.ndarray | None = None  # scaled

    @property
    def n(self) -> int:
        return len(self.t)

    def target(self) -> np.ndarray:
        return self.t >= self.tstart


def _floats(line: str) -> list[float]:
    return [float(v.replace("D", "E").replace("d", "e")) for v in line.split()]


def _aniso(dd, ee, ff):
    s = np.sqrt(1 - ff ** 2)
    return ee / dd / s, ff / s, dd / ee / s


def read_conf(path: str) -> dict:
    """hist-etas5pa.conf / hist-etas-mk.conf: lines 2 and 3 only."""
    with open(path) as fh:
        lines = fh.readlines()
    tx, ty, tz, nn = _floats(lines[1].split("!")[0])[:4]
    xmin, ymin, xmg0, zmin, tstart, bi2 = _floats(lines[2].split("!")[0])[:6]
    return dict(tx=tx, ty=ty, tz=tz, nn=int(nn), xmin=xmin, ymin=ymin,
                xmg0=xmg0, zmin=zmin, tstart=tstart, bi2=bi2)


def read_delone2(path: str, conf: dict) -> Catalog:
    with open(path) as fh:
        head = fh.readline().split()
        npnt, npext, ndelo = int(head[0]), int(head[1]), int(head[2])
        rows = np.array([_floats(fh.readline())[1:8] for _ in range(npext)])
        trows = np.array([_floats(fh.readline()) for _ in range(ndelo)])
    cf = np.cos((conf["ymin"] + conf["ty"] / 2) * FORTRAN_PI / 180.0)
    ev = rows[:npnt]
    qa, qb, qc = _aniso(ev[:, 4], ev[:, 5], ev[:, 6])
    return Catalog(
        x=ev[:, 0] * cf, y=ev[:, 1].copy(), mag=ev[:, 2] - conf["xmg0"],
        t=ev[:, 3] - conf["zmin"], qa=qa, qb=qb, qc=qc,
        tx=conf["tx"] * cf, ty=conf["ty"], tz=conf["tz"], tstart=conf["tstart"],
        xmg0=conf["xmg0"], xmin=conf["xmin"], ymin=conf["ymin"], cosfac=cf,
        node_x=rows[:, 0] * cf, node_y=rows[:, 1].copy(),
        tri=trows[:, 1:4].astype(int) - 1, tri_area=trows[:, 4] * cf)


def read_out3(path: str, conf: dict) -> Catalog:
    """etas2aniso.out3, read the way st-etas reads it: absolute lon/lat."""
    with open(path) as fh:
        fh.readline()
        rows = np.array([_floats(l)[1:8] for l in fh if l.strip()])
    rows = rows[(rows[:, 2] >= conf["xmg0"]) & (rows[:, 1] != conf["ymin"])]
    cf = np.cos((conf["ymin"] + conf["ty"] / 2) * FORTRAN_PI / 180.0)
    qa, qb, qc = _aniso(rows[:, 4], rows[:, 5], rows[:, 6])
    return Catalog(
        x=rows[:, 0] * cf, y=rows[:, 1].copy(), mag=rows[:, 2] - conf["xmg0"],
        t=rows[:, 3] - conf["zmin"], qa=qa, qb=qb, qc=qc,
        tx=conf["tx"] * cf, ty=conf["ty"], tz=conf["tz"] - conf["zmin"],
        tstart=conf["tstart"], xmg0=conf["xmg0"], xmin=conf["xmin"],
        ymin=conf["ymin"], cosfac=cf)


# ---------------------------------------------------------------- parameters
@dataclass
class Fields:
    """Parameter values at the Delaunay nodes (events first).

    ``a`` is the halved alpha used inside the kernel; ``alpha = 2a``.
    """
    mu: np.ndarray
    K: np.ndarray
    a: np.ndarray
    p: np.ndarray
    q: np.ndarray
    c: float
    d: float
    const_mu: bool = False  # background integral = mu * area (st-etas)

    def at_events(self, n: int) -> "Fields":
        return replace(self, mu=self.mu[:n], K=self.K[:n], a=self.a[:n],
                       p=self.p[:n], q=self.q[:n])


def _upda_tokens(path: str) -> list[float]:
    with open(path) as fh:
        return _floats(fh.read())


def read_5pa(path: str) -> Fields:
    v = _upda_tokens(path)
    m = int(v[6])
    a1, a2, a3, aa4, a5, a6, a7 = v[8:15]
    xd = np.array(v[15:15 + m]).reshape(5, m // 5)
    return Fields(mu=a1 * np.exp(xd[0]), K=a2 * np.exp(xd[1]),
                  a=aa4 / 2 * np.exp(xd[2]), p=a5 * np.exp(xd[3]),
                  q=a7 * np.exp(xd[4]), c=a3, d=a6)


def read_mk(path: str) -> Fields:
    v = _upda_tokens(path)
    m = int(v[3])
    a1, a2, a3, aa4, a5, a6, a7 = v[4:11]
    xd = np.array(v[11:11 + m]).reshape(2, m // 2)
    one = np.ones(m // 2)
    return Fields(mu=a1 * np.exp(xd[0]), K=a2 * np.exp(xd[1]),
                  a=aa4 / 2 * one, p=a5 * one, q=a7 * one, c=a3, d=a6)


def constant_fields(theta, n: int) -> Fields:
    """theta = (mu, K, c, alpha, p, d, q) in the conf/upda order (unhalved alpha)."""
    mu, K, c, alpha, p, d, q = theta
    one = np.ones(n)
    return Fields(mu=mu * one, K=K * one, a=alpha / 2 * one, p=p * one,
                  q=q * one, c=c, d=d, const_mu=True)


# ---------------------------------------------------------------- triggering
@dataclass
class Pairs:
    i: np.ndarray   # triggered
    j: np.ndarray   # trigger, j < i
    dt: np.ndarray
    dx: np.ndarray
    dy: np.ndarray


def neighbour_pairs(cat: Catalog, bi2: float = 2.0, rule: str = "circle") -> Pairs:
    """Candidate (i, j) pairs, j < i, inside the programs' search radius.

    rmax_j = (0.015 * 10^(0.5 M_j - 2) + 0.3) * bi2 degrees; hist-etas uses a
    circle.  rule='all' keeps every earlier event, which is what st-etas's
    func17 does (it builds a square neighbour list but never uses it).
    """
    n = cat.n
    dx = cat.x[:, None] - cat.x[None, :]
    dy = cat.y[:, None] - cat.y[None, :]
    rmax = (0.015 * 10.0 ** (0.5 * (cat.mag + cat.xmg0) - 2.0) + 0.3) * bi2
    if rule == "circle":
        keep = dx ** 2 + dy ** 2 <= rmax[None, :] ** 2
    else:
        keep = np.ones((n, n), dtype=bool)
    keep &= np.tri(n, n, -1, dtype=bool)
    i, j = np.nonzero(keep)
    return Pairs(i=i, j=j, dt=cat.t[i] - cat.t[j], dx=dx[i, j], dy=dy[i, j])


def trig_rate(cat: Catalog, f: Fields, pr: Pairs) -> np.ndarray:
    """Per-pair triggering contribution to lambda at event i."""
    j = pr.j
    s = np.exp(-f.a[j] * cat.mag[j])
    u, v = pr.dx * s, pr.dy * s
    Q = cat.qa[j] * u ** 2 - 2 * cat.qb[j] * u * v + cat.qc[j] * v ** 2
    return f.K[j] / (pr.dt + f.c) ** f.p[j] * (Q + f.d) ** (-f.q[j])


def intensity_at_events(cat: Catalog, f: Fields, pr: Pairs) -> tuple[np.ndarray, np.ndarray]:
    """(background, triggered) parts of lambda at every event."""
    trig = np.bincount(pr.i, weights=trig_rate(cat, f, pr), minlength=cat.n)
    return f.mu[:cat.n].copy(), trig


# ---------------------------------------------------------------- compensator
def omori_int(u, c, p):
    """int_0^u (s + c)^-p ds."""
    u = np.maximum(u, 0.0)
    return ((u + c) ** (1 - p) - c ** (1 - p)) / (1 - p)


def region_mass(cat: Catalog, f: Fields, nseg: int = 128,
                fortran_offsets: bool = False, clamp: bool = True) -> np.ndarray:
    """Spatial integral of each event's kernel over the box [0,tx] x [0,ty].

    The box is cut into rays from the event to nseg segments per side; each
    sector is replaced by a disc sector whose squared Q-radius is the mean of
    the two end values, as in the programs.  ``fortran_offsets`` keeps the
    programs' (+x_i, +y_i) offsets on the left and bottom sides, which flip
    the sign of the cross term of Q there.
    """
    x, y, qa, qb, qc = cat.x, cat.y, cat.qa, cat.qb, cat.qc
    tx, ty = cat.tx, cat.ty
    bx = np.linspace(0.0, tx, nseg + 1)
    by = np.linspace(0.0, ty, nseg + 1)
    a = f.a[:cat.n]
    q = f.q[:cat.n]
    scl = np.exp(2 * a * cat.mag)
    d = f.d

    def Q(u, v):
        return qa[:, None] * u ** 2 - 2 * qb[:, None] * u * v + qc[:, None] * v ** 2

    sides = []
    xl = x[:, None] if fortran_offsets else -x[:, None]
    yl = y[:, None] if fortran_offsets else -y[:, None]
    # left, top, right, bottom -- (u1, v1, u2, v2, chord^2, skip-mask)
    sides.append((xl + 0 * by[None, :-1], by[None, :-1] - y[:, None],
                  xl + 0 * by[None, 1:], by[None, 1:] - y[:, None],
                  (ty / nseg) ** 2, x <= 0.0))
    sides.append((bx[None, :-1] - x[:, None], ty - y[:, None] + 0 * bx[None, :-1],
                  bx[None, 1:] - x[:, None], ty - y[:, None] + 0 * bx[None, 1:],
                  (tx / nseg) ** 2, y >= ty))
    sides.append((tx - x[:, None] + 0 * by[None, :-1], by[None, :-1] - y[:, None],
                  tx - x[:, None] + 0 * by[None, 1:], by[None, 1:] - y[:, None],
                  (ty / nseg) ** 2, x >= tx))
    sides.append((bx[None, :-1] - x[:, None], yl + 0 * bx[None, :-1],
                  bx[None, 1:] - x[:, None], yl + 0 * bx[None, 1:],
                  (tx / nseg) ** 2, y <= 0.0))

    mass = np.zeros(cat.n)
    for u1, v1, u2, v2, rr3, skip in sides:
        rrd = (Q(u1, v1) + Q(u2, v2)) / 2
        rr1 = u1 ** 2 + v1 ** 2
        rr2 = u2 ** 2 + v2 ** 2
        cosd = (rr1 + rr2 - rr3) / (2 * np.sqrt(rr1) * np.sqrt(rr2))
        if clamp:
            cosd = np.clip(cosd, -1.0, 1.0)
        delta = np.arccos(cosd)
        sp = np.pi / (1 - q[:, None]) * scl[:, None] * (
            (rrd / scl[:, None] + d) ** (1 - q[:, None]) - d ** (1 - q[:, None]))
        contrib = (delta / (2 * np.pi) * sp).sum(axis=1)
        mass += np.where(skip, 0.0, contrib)
    return mass


def background_integral(cat: Catalog, f: Fields) -> float:
    """int_A mu(x, y) dx dy, per day."""
    if f.const_mu:
        return float(f.mu[0]) * cat.tx * cat.ty
    tri = cat.tri
    return float(np.sum(cat.tri_area * f.mu[tri].sum(axis=1) / 3))


def _st_volume(cat: Catalog, nox: int = 101, noy: int = 101) -> float:
    """st-etas's trapezoid 'volume' of exp(0) with its corner weights."""
    w = np.ones((nox, noy))
    w[[0, -1], :] = 0.5
    w[:, [0, -1]] = 0.5
    w[0, 0] = w[-1, -1] = 1 / 3
    w[0, -1] = w[-1, 0] = 1 / 6
    return w.sum() * (cat.tx / (nox - 1)) * (cat.ty / (noy - 1))


# ---------------------------------------------------------------- log-likelihood
def loglik(cat: Catalog, f: Fields, pr: Pairs, mode: str = "exact",
           program: str = "hist", nseg: int | None = None, mass=None) -> dict:
    """Log-likelihood over the target period.

    program = 'hist' (hist-etas5pa / hist-etas-mk) or 'st' (st-etas), only
    relevant in fortran mode.
    """
    fe = f.at_events(cat.n)
    bg, trig = intensity_at_events(cat, fe, pr)
    lam = bg + trig
    tgt = cat.target()
    if mode == "fortran":
        tgt = tgt & (np.arange(cat.n) >= 1)
        # every program adds log(mu) of the first event, precursory or not
        sum_log = np.log(bg[0]) + np.sum(np.log(lam[tgt]))
        if program == "st":
            S = region_mass(cat, fe, nseg=nseg or 128, fortran_offsets=True, clamp=False)
            late = cat.t >= cat.tstart
            trig_int = np.sum((fe.K * S * omori_int(cat.tz - cat.t, fe.c, fe.p))[late])
            bg_int = fe.mu[0] * cat.tz * _st_volume(cat)
        else:
            S = region_mass(cat, fe, nseg=nseg or 8, fortran_offsets=True)
            T = omori_int(cat.tz - cat.t, fe.c, fe.p) - omori_int(cat.tstart - cat.t, fe.c, fe.p)
            trig_int = np.sum(fe.K * S * T)
            bg_int = background_integral(cat, f) * (cat.tz - cat.tstart)
    else:
        sum_log = np.sum(np.log(lam[tgt]))
        S = region_mass(cat, fe, nseg=nseg or 128) if mass is None else mass
        T = omori_int(cat.tz - cat.t, fe.c, fe.p) - omori_int(cat.tstart - cat.t, fe.c, fe.p)
        trig_int = np.sum(fe.K * S * T)
        bg_int = background_integral(cat, f) * (cat.tz - cat.tstart)
    return dict(loglik=sum_log - trig_int - bg_int, sum_log=sum_log,
                trig_int=trig_int, bg_int=bg_int, n_target=int(tgt.sum()))


# ---------------------------------------------------------------- compensator path
def compensator_at_events(cat: Catalog, f: Fields, mass: np.ndarray) -> np.ndarray:
    """Lambda(t_i) = int_tstart^t_i int_A lambda, at every target event (exact mode)."""
    fe = f.at_events(cat.n)
    tgt = np.nonzero(cat.target())[0]
    ti = cat.t[tgt]
    B = background_integral(cat, f)
    # triggering: sum over all earlier events j of K_j S_j [G(t_i - t_j) - G(tstart - t_j)]
    w = fe.K * mass
    g0 = omori_int(cat.tstart - cat.t, fe.c, fe.p)
    lam = B * (ti - cat.tstart)
    for k, (idx, t) in enumerate(zip(tgt, ti)):
        j = slice(0, idx)
        lam[k] += np.sum(w[j] * (omori_int(t - cat.t[j], fe.c, fe.p[j]) - g0[j]))
    return lam
