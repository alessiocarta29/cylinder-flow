#!/usr/bin/env python3
# Ricalcola la forza viscosa sul cilindro con una formula a tre punti per il
# gradiente di velocita' a parete (parete + prime due celle), invece della
# differenza a due punti usata da OpenFOAM.
#
# Uso:  python3 scripts/wall_shear.py runs/re20/L1
# Prima, nella cartella del caso:  postProcess -func writeCellCentres -latestTime

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wall_pressure import read_entry, read_field, latest_time, wall_faces

def read_nu(path):
    # la riga e':  nu  [0 2 -1 0 0 0 0] 0.05;
    # read_entry non va bene qui, perche' prende il secondo campo della riga
    for line in open(path):
        parts = line.split()
        if parts and parts[0] == "nu":
            for token in parts:
                if token.endswith(";"):
                    return float(token.rstrip(";"))
    sys.exit(f"{path}: voce 'nu' non trovata")

def drag_coefficients(case):
    n_r = int(read_entry(os.path.join(case, "system/blockMeshDict"), "nR"))
    n_theta = int(read_entry(os.path.join(case, "system/blockMeshDict"), "nTheta"))
    nu = read_nu(os.path.join(case, "constant/transportProperties"))
    time_dir = latest_time(case)
    u = read_field(os.path.join(time_dir, "U"))
    centres = read_field(os.path.join(time_dir, "C"))
    if len(u) != 4 * n_r * n_theta or len(centres) != len(u):
        sys.exit(f"{time_dir}: attese {4 * n_r * n_theta} celle, trovate {len(u)}")

    faces = wall_faces(n_theta, 0.5, 1.0)
    fx_two = 0.0      # gradiente da due punti (parete e prima cella)
    fx_three = 0.0    # gradiente da tre punti (parete e prime due celle)

    for k, (centre, normal, area) in enumerate(faces):
        block, j = divmod(k, n_theta)
        cell1 = block * n_r * n_theta + j * n_r
        cell2 = cell1 + 1

        # tangente alla parete, ruotata di 90 gradi rispetto alla normale
        tangent = (-normal[1], normal[0])

        # distanze dei due centri cella dalla faccia, lungo la normale
        d1 = ((centres[cell1][0] - centre[0]) * normal[0]
              + (centres[cell1][1] - centre[1]) * normal[1])
        d2 = ((centres[cell2][0] - centre[0]) * normal[0]
              + (centres[cell2][1] - centre[1]) * normal[1])

        # componenti tangenziali della velocita' nelle due celle
        ut1 = u[cell1][0] * tangent[0] + u[cell1][1] * tangent[1]
        ut2 = u[cell2][0] * tangent[0] + u[cell2][1] * tangent[1]

        # a parete la velocita' e' nulla (aderenza), quindi il terzo punto vale 0
        grad_two = ut1 / d1
        grad_three = (ut1 * d2 * d2 - ut2 * d1 * d1) / (d1 * d2 * (d2 - d1))

        # forza sul cilindro: F = somma di nu * (du_t/dn) * tangente * area
        fx_two += nu * grad_two * tangent[0] * area
        fx_three += nu * grad_three * tangent[0] * area

    # rho = 1, U = 1, A = D * span = 1  ->  C = 2 F  (eq. 9.98)
    return 2.0 * fx_two, 2.0 * fx_three

if __name__ == "__main__":
    case = sys.argv[1]
    cd_two, cd_three = drag_coefficients(case)
    print(f"{case}")
    print(f"  Cd viscosa, due punti          : {cd_two:.8f}")
    print(f"  Cd viscosa, tre punti          : {cd_three:.8f}")
    print(f"  differenza                     : {cd_three - cd_two:+.8f}")
