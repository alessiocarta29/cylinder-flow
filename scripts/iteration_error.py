#!/usr/bin/env python3
"""Iteration-error estimate on Cd for a steady SIMPLE run.

Ferziger, Peric & Street, sec. 5.7:
    delta_n = Cd_{n+1} - Cd_n
    lambda  ~ |delta_n| / |delta_{n-1}|     (averaged over the last iterations)
    eps     ~ |delta_n| / (1 - lambda)      (valid for a real dominant eigenvalue)

Usage (from the case directory):
    python3 ../scripts/iteration_error.py [path/to/coefficient.dat]
"""
import sys

path = sys.argv[1] if len(sys.argv) > 1 else "postProcessing/forceCoeffs/0/coefficient.dat"
window = 10  # number of iterations used to average lambda

it, cd = [], []
with open(path) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split()
        it.append(int(float(cols[0])))
        cd.append(float(cols[1]))

delta = [b - a for a, b in zip(cd[:-1], cd[1:])]
tail = delta[-(window + 1):]

ratios = [abs(d1) / abs(d0) for d0, d1 in zip(tail[:-1], tail[1:]) if d0 != 0.0]
sign_changes = sum(1 for d0, d1 in zip(tail[:-1], tail[1:]) if d0 * d1 < 0.0)

print(f"iterations           : {it[-1]}")
print(f"Cd (last)            : {cd[-1]:.12f}")
print("last deltas          : " + " ".join(f"{d:+.2e}" for d in delta[-6:]))
print(f"sign changes (last {window}): {sign_changes}")

if not ratios or delta[-1] == 0.0:
    print("lambda               : not computable (last deltas below write resolution:")
    print("                       the iteration error is below ~1e-12 or writePrecision is too low)")
    sys.exit(0)

lam = sum(ratios) / len(ratios)
print(f"lambda (mean)        : {lam:.4f}")

if sign_changes > window // 2:
    print("iteration error      : not estimated (oscillating deltas, complex eigenvalues likely;")
    print("                       the real-eigenvalue formula does not apply, see sec. 5.7)")
    print(f"                       max |delta| over last {window}: {max(abs(d) for d in tail):.2e}")
elif lam >= 1.0:
    print("iteration error      : not estimated (lambda >= 1: not converging)")
else:
    eps = abs(delta[-1]) / (1.0 - lam)
    print(f"iteration error Cd   : {eps:.2e}  ({eps / abs(cd[-1]):.1e} relative)")
