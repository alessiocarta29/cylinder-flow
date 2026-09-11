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

DEFAULT_PATH = "postProcessing/forceCoeffs/0/coefficient.dat"


def read_cd(path):
    """Iteration numbers and Cd values from a forceCoeffs coefficient.dat file."""
    it, cd = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            it.append(int(float(cols[0])))
            cd.append(float(cols[1]))
    return it, cd


def estimate(path, window=10):
    """Sec. 5.7 estimate. 'status' is 'ok' only when the real-eigenvalue formula applies."""
    it, cd = read_cd(path)
    result = {
        "iterations": it[-1] if it else 0,
        "cd": cd[-1] if cd else None,
        "last_deltas": [],
        "sign_changes": 0,
        "window": window,
        "lambda": None,
        "eps": None,
        "eps_rel": None,
        "max_tail_delta": None,
    }
    if len(cd) < 3:
        result["status"] = "too_few_iterations"
        return result

    delta = [b - a for a, b in zip(cd[:-1], cd[1:])]
    tail = delta[-(window + 1):]
    ratios = [abs(d1) / abs(d0) for d0, d1 in zip(tail[:-1], tail[1:]) if d0 != 0.0]
    result["last_deltas"] = delta[-6:]
    result["sign_changes"] = sum(1 for d0, d1 in zip(tail[:-1], tail[1:]) if d0 * d1 < 0.0)

    if not ratios or delta[-1] == 0.0:
        result["status"] = "below_resolution"
        return result

    lam = sum(ratios) / len(ratios)
    result["lambda"] = lam

    if result["sign_changes"] > window // 2:
        result["status"] = "oscillating"
        result["max_tail_delta"] = max(abs(d) for d in tail)
    elif lam >= 1.0:
        result["status"] = "not_converging"
    else:
        eps = abs(delta[-1]) / (1.0 - lam)
        result.update(status="ok", eps=eps, eps_rel=eps / abs(cd[-1]))
    return result


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_PATH
    r = estimate(path)

    print(f"iterations           : {r['iterations']}")
    if r["status"] == "too_few_iterations":
        print("iteration error      : not estimated (fewer than 3 iterations)")
        return
    print(f"Cd (last)            : {r['cd']:.12f}")
    print("last deltas          : " + " ".join(f"{d:+.2e}" for d in r["last_deltas"]))
    print(f"sign changes (last {r['window']}): {r['sign_changes']}")

    if r["status"] == "below_resolution":
        print("lambda               : not computable (last deltas below write resolution:")
        print("                       the iteration error is below ~1e-12 or writePrecision is too low)")
        return

    print(f"lambda (mean)        : {r['lambda']:.4f}")
    if r["status"] == "oscillating":
        print("iteration error      : not estimated (oscillating deltas, complex eigenvalues likely;")
        print("                       the real-eigenvalue formula does not apply, see sec. 5.7)")
        print(f"                       max |delta| over last {r['window']}: {r['max_tail_delta']:.2e}")
    elif r["status"] == "not_converging":
        print("iteration error      : not estimated (lambda >= 1: not converging)")
    else:
        print(f"iteration error Cd   : {r['eps']:.2e}  ({r['eps_rel']:.1e} relative)")


if __name__ == "__main__":
    main()
