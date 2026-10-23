from __future__ import annotations
import math
from fractions import Fraction
import numpy as np


def symplectic_J(n: int = 3) -> np.ndarray:
    I = np.eye(n, dtype=int)
    Z = np.zeros((n, n), dtype=int)
    return np.block([[Z, I], [-I, Z]])


def pairing(x: np.ndarray, y: np.ndarray, J: np.ndarray) -> int:
    return int(np.asarray(x, dtype=object).T @ np.asarray(J, dtype=object) @ np.asarray(y, dtype=object))


def pairings_to_basis(x: np.ndarray, basis: np.ndarray, J: np.ndarray) -> tuple[int, ...]:
    # basis columns are reference charge vectors.
    return tuple(pairing(x, basis[:, i], J) for i in range(basis.shape[1]))


def random_unimodular(n: int, rng: np.random.Generator, steps: int = 8) -> np.ndarray:
    A = np.eye(n, dtype=int)
    for _ in range(steps):
        op = int(rng.integers(0, 3))
        if op == 0:
            i, j = rng.choice(n, size=2, replace=False)
            A[[i, j], :] = A[[j, i], :]
        elif op == 1:
            i = int(rng.integers(n))
            A[i, :] *= -1
        else:
            i, j = rng.choice(n, size=2, replace=False)
            k = int(rng.integers(-3, 4))
            if k == 0:
                k = 1
            A[i, :] += k * A[j, :]
    det = round(np.linalg.det(A))
    if abs(det) != 1:
        raise RuntimeError(f"Generated non-unimodular matrix det={det}")
    return A


def random_symmetric_int(n: int, rng: np.random.Generator, bound: int = 2) -> np.ndarray:
    M = rng.integers(-bound, bound + 1, size=(n, n))
    return np.triu(M) + np.triu(M, 1).T


def random_symplectic(n: int, rng: np.random.Generator, stages: int = 3) -> np.ndarray:
    I = np.eye(n, dtype=int)
    Z = np.zeros((n, n), dtype=int)
    U = np.eye(2 * n, dtype=int)
    for _ in range(stages):
        A = random_unimodular(n, rng, steps=5)
        AinvT = np.rint(np.linalg.inv(A).T).astype(int)
        UA = np.block([[A, Z], [Z, AinvT]])
        B = random_symmetric_int(n, rng)
        UB = np.block([[I, B], [Z, I]])
        C = random_symmetric_int(n, rng)
        UC = np.block([[I, Z], [C, I]])
        U = UC @ UB @ UA @ U
    J = symplectic_J(n)
    if not np.array_equal(U.T @ J @ U, J):
        raise RuntimeError("Symplectic construction failed")
    return U


def quintic_euler_line_bundle(d: int) -> Fraction:
    # On a smooth quintic CY3: H^3=5, c2.H=50,
    # chi(O_X(d)) = (5 d^3 + 25 d)/6.
    return Fraction(5 * d**3 + 25 * d, 6)


def hrr_pairing_line_bundles(a: int, b: int) -> Fraction:
    return quintic_euler_line_bundle(b - a)


def random_unitary(n: int, rng: np.random.Generator) -> np.ndarray:
    Z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Q, R = np.linalg.qr(Z)
    phase = np.diag(R)
    phase = np.where(np.abs(phase) == 0, 1.0, phase / np.abs(phase))
    Q = Q @ np.diag(np.conjugate(phase))
    return Q


def random_su3_hermitian(rng: np.random.Generator) -> np.ndarray:
    Z = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    H = (Z + Z.conjugate().T) / 2.0
    H -= np.trace(H) / 3.0 * np.eye(3)
    return H


def gauge_invariants(H: np.ndarray) -> np.ndarray:
    eig = np.linalg.eigvalsh(H)
    c2 = np.trace(H @ H).real
    c3 = np.trace(H @ H @ H).real
    return np.array([c2, c3, *eig], dtype=float)


def wilson_interval(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return (float('nan'), float('nan'))
    p = k / n
    den = 1 + z*z/n
    center = (p + z*z/(2*n)) / den
    half = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return center - half, center + half
