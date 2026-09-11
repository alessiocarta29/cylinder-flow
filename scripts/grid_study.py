#!/usr/bin/env python3
"""Grid-refinement study for the 2D cylinder at Re = 20.

Builds the five O-grids of Ferziger, Peric & Street (sec. 9.12) from the
template case re20/, runs blockMesh, checkMesh and simpleFoam on each level,
and collects forces and iteration errors in results/re20_grid_study.csv.

Only the grid changes between levels. For level k = 1..5:
    cells per block around the cylinder   nTheta = 6 * 2^(k-1)   (4 blocks)
    radial cells                          nR     = 16 * 2^(k-1)
    radial cell-to-cell ratio             r_k    = 1.25^(1 / 2^(k-1))
    blockMesh radial grading              r_k^(nR - 1)   (last/first cell)
Iteration cap and write interval are run controls, not numerics.

Usage (OpenFOAM environment loaded):
    python3 scripts/grid_study.py          # all five levels
    python3 scripts/grid_study.py 1 2      # selected levels only
Existing rows of other levels in the CSV are kept.
"""
import argparse
import csv
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

from iteration_error import estimate

REPO = Path(__file__).resolve().parents[1]
TEMPLATE = REPO / "re20"
RUNS = REPO / "runs" / "re20"
RESULTS = REPO / "results" / "re20_grid_study.csv"

N_THETA_1 = 6          # level 1: cells per block around the cylinder
N_R_1 = 16             # level 1: radial cells
R_1 = 1.25             # level 1: radial cell-to-cell ratio

END_TIME_CAP = 300000  # runs normally stop earlier on residualControl
WRITE_INTERVAL = 10000

# forceCoeffs uses U = 1, Aref = D * span = 1 and rho = 1, so C = 2 F (eq. 9.98)
FORCE_TO_COEFF = 2.0

FIELDS = [
    "level", "n_theta_total", "n_r", "cells", "radial_ratio", "radial_grading",
    "max_non_orthogonality", "iterations", "converged",
    "cd", "cd_pressure", "cd_viscous", "cl",
    "iteration_error_rel", "iteration_error_status", "wall_time_s", "commit",
]


def level_parameters(k):
    n_theta = N_THETA_1 * 2 ** (k - 1)
    n_r = N_R_1 * 2 ** (k - 1)
    r = R_1 ** (1.0 / 2 ** (k - 1))
    return n_theta, n_r, r, r ** (n_r - 1)


def set_entry(path, keyword, value):
    """Replace a top-level 'keyword value;' entry; exactly one match is required."""
    text = path.read_text()
    pattern = rf"^({re.escape(keyword)}\s+)[^;\s]+;"
    new, count = re.subn(pattern, rf"\g<1>{value};", text, flags=re.MULTILINE)
    if count != 1:
        sys.exit(f"{path}: expected one top-level '{keyword}' entry, found {count}")
    path.write_text(new)


def prepare_case(k):
    case = RUNS / f"L{k}"
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    for sub in ("0", "constant", "system"):
        shutil.copytree(TEMPLATE / sub, case / sub, ignore=shutil.ignore_patterns("polyMesh"))
    return case


def run(case, app):
    log = case / f"log.{app}"
    with open(log, "w") as f:
        proc = subprocess.run([app], cwd=case, stdout=f, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        sys.exit(f"{app} failed for {case.name}: see {log}")
    return log


def file_contains(path, text):
    with open(path) as f:
        return any(text in line for line in f)


def last_data_row(path):
    last = None
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                last = line
    if last is None:
        sys.exit(f"{path}: no data")
    return [float(x) for x in last.split()]


def git_commit():
    try:
        out = subprocess.run(["git", "describe", "--always", "--dirty"], cwd=REPO,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def read_rows():
    if not RESULTS.exists():
        return {}
    with open(RESULTS, newline="") as f:
        return {int(row["level"]): row for row in csv.DictReader(f)}


def write_rows(rows):
    RESULTS.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for k in sorted(rows):
            writer.writerow(rows[k])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("levels", nargs="*", type=int, help="grid levels (1-5), default all")
    levels = parser.parse_args().levels or [1, 2, 3, 4, 5]
    if any(k not in range(1, 6) for k in levels):
        sys.exit("levels must be between 1 and 5")

    for app in ("blockMesh", "checkMesh", "simpleFoam"):
        if shutil.which(app) is None:
            sys.exit(f"{app} not found: load the OpenFOAM environment first")

    commit = git_commit()
    if commit.endswith("-dirty"):
        print(f"Warning: uncommitted changes in the repository; results tagged {commit}\n")

    rows = read_rows()

    for k in levels:
        n_theta, n_r, r, grading = level_parameters(k)
        case = prepare_case(k)

        bmd = case / "system" / "blockMeshDict"
        set_entry(bmd, "nTheta", n_theta)
        set_entry(bmd, "nR", n_r)
        set_entry(bmd, "radialGrading", f"{grading:.12g}")
        cdict = case / "system" / "controlDict"
        set_entry(cdict, "endTime", END_TIME_CAP)
        set_entry(cdict, "writeInterval", WRITE_INTERVAL)

        print(f"L{k}: {4 * n_theta}x{n_r}, {4 * n_theta * n_r} cells, running...", flush=True)
        start = time.time()
        run(case, "blockMesh")
        check_log = run(case, "checkMesh")
        if not file_contains(check_log, "Mesh OK."):
            sys.exit(f"checkMesh did not report 'Mesh OK.' for L{k}: see {check_log}")
        match = re.search(r"Mesh non-orthogonality Max:\s*([-+0-9.eE]+)", check_log.read_text())
        solver_log = run(case, "simpleFoam")
        wall = time.time() - start

        pp = case / "postProcessing"
        force = last_data_row(pp / "forces" / "0" / "force.dat")
        coeff = last_data_row(pp / "forceCoeffs" / "0" / "coefficient.dat")
        cd, cl = coeff[1], coeff[2]
        if abs(cd - FORCE_TO_COEFF * force[1]) > 1e-9 * abs(cd):
            print(f"  Warning: Cd from forceCoeffs ({cd}) differs from 2*Fx ({FORCE_TO_COEFF * force[1]})")

        ie = estimate(pp / "forceCoeffs" / "0" / "coefficient.dat")
        converged = file_contains(solver_log, "SIMPLE solution converged")

        rows[k] = {
            "level": k,
            "n_theta_total": 4 * n_theta,
            "n_r": n_r,
            "cells": 4 * n_theta * n_r,
            "radial_ratio": f"{r:.12g}",
            "radial_grading": f"{grading:.12g}",
            "max_non_orthogonality": match.group(1) if match else "",
            "iterations": int(force[0]),
            "converged": "yes" if converged else "no",
            "cd": f"{cd:.12g}",
            "cd_pressure": f"{FORCE_TO_COEFF * force[4]:.12g}",
            "cd_viscous": f"{FORCE_TO_COEFF * force[7]:.12g}",
            "cl": f"{cl:.3e}",
            "iteration_error_rel": f"{ie['eps_rel']:.2e}" if ie["status"] == "ok" else "",
            "iteration_error_status": ie["status"],
            "wall_time_s": f"{wall:.1f}",
            "commit": commit,
        }
        write_rows(rows)

        err = rows[k]["iteration_error_rel"] or ie["status"]
        print(f"  iterations {int(force[0])} (converged: {rows[k]['converged']}), "
              f"Cd {cd:.8f} = {FORCE_TO_COEFF * force[4]:.8f} (p) + {FORCE_TO_COEFF * force[7]:.8f} (v), "
              f"iteration error {err}, {wall:.0f} s", flush=True)

    print(f"\nResults: {RESULTS.relative_to(REPO)}")


if __name__ == "__main__":
    main()
