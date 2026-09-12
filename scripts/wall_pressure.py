#!/usr/bin/env python3
# Ricalcola la forza di pressione sul cilindro estrapolando linearmente
# la pressione di parete dalle prime due celle radiali (Ferziger, sez. 7.1),
# invece di usare il valore della prima cella (zeroGradient).
#
# Uso:  python3 scripts/wall_pressure.py runs/re20/L1
# Prima, nella cartella del caso:  postProcess -func writeCellCentres -latestTime

import math
import os
import sys

# ---------- lettura dei file OpenFOAM ----------

def read_entry(path, keyword):
    # legge una voce "keyword valore;" dal file
    for line in open(path):
        parts = line.split()
        if len(parts) >= 2 and parts[0] == keyword:
            return float(parts[1].rstrip(";"))
    sys.exit(f"{path}: voce '{keyword}' non trovata")

def read_field(path):
    # legge internalField di un campo scalare o vettoriale;
    # restituisce una lista di numeri o di terne
    text = open(path).read()
    start = text.index("internalField")
    is_vector = "List<vector>" in text[start:start + 100]
    open_par = text.index("(", start)              # parentesi che apre la lista
    close_par = text.index("\n)", open_par)        # ")" a inizio riga che la chiude
    body = text[open_par + 1:close_par]

    if not is_vector:
        return [float(t) for t in body.split()]

    values = []
    for line in body.split("\n"):
        line = line.strip()
        if line.startswith("("):
            x, y, z = line.strip("()").split()
            values.append((float(x), float(y), float(z)))
    return values

def latest_time(case):
    # cartella di tempo piu' alta, escludendo 0
    times = []
    for name in os.listdir(case):
        try:
            t = float(name)
        except ValueError:
            continue
        if t > 0 and os.path.isdir(os.path.join(case, name)):
            times.append((t, name))
    if not times:
        sys.exit(f"{case}: nessuna cartella di tempo")
    return os.path.join(case, max(times)[1])

# ---------- geometria della mesh a O ----------
# 4 blocchi, ciascuno con nR celle radiali (indice i, la 0 tocca il cilindro)
# e nTheta celle circonferenziali (indice j). blockMesh numera le celle con i
# piu' veloce, poi j, poi blocco: cella = blocco*nR*nTheta + j*nR + i.

def wall_faces(n_theta, radius, span):
    # per ogni faccia sul cilindro: centro, normale uscente dal cilindro, area
    faces = []
    d_theta = 0.5 * math.pi / n_theta          # 90 gradi divisi per blocco
    for block in range(4):
        theta0 = 1.25 * math.pi + block * 0.5 * math.pi   # si parte da 225 gradi
        for j in range(n_theta):
            theta = theta0 + (j + 0.5) * d_theta          # angolo di meta' faccia
            normal = (math.cos(theta), math.sin(theta))
            chord = 2.0 * radius * math.sin(0.5 * d_theta)
            centre = (radius * math.cos(theta) * math.cos(0.5 * d_theta),
                      radius * math.sin(theta) * math.cos(0.5 * d_theta))
            faces.append((centre, normal, chord * span))
    return faces

# ---------- calcolo ----------

def drag_coefficients(case):
    n_r = int(read_entry(os.path.join(case, "system/blockMeshDict"), "nR"))
    n_theta = int(read_entry(os.path.join(case, "system/blockMeshDict"), "nTheta"))
    time_dir = latest_time(case)
    p = read_field(os.path.join(time_dir, "p"))
    centres = read_field(os.path.join(time_dir, "C"))
    if len(p) != 4 * n_r * n_theta or len(centres) != len(p):
        sys.exit(f"{time_dir}: attese {4 * n_r * n_theta} celle, trovate {len(p)} e {len(centres)}")

    faces = wall_faces(n_theta, 0.5, 1.0)
    fx_first = 0.0      # pressione di parete = valore della prima cella
    fx_linear = 0.0     # pressione di parete estrapolata dalle prime due

    for k, (centre, normal, area) in enumerate(faces):
        block, j = divmod(k, n_theta)
        cell1 = block * n_r * n_theta + j * n_r
        cell2 = cell1 + 1

        # distanze dei due centri cella dalla faccia, lungo la normale
        d1 = ((centres[cell1][0] - centre[0]) * normal[0]
              + (centres[cell1][1] - centre[1]) * normal[1])
        d2 = ((centres[cell2][0] - centre[0]) * normal[0]
              + (centres[cell2][1] - centre[1]) * normal[1])

        p_linear = p[cell1] - d1 * (p[cell2] - p[cell1]) / (d2 - d1)

        # forza sul cilindro: F = - somma di p * normale * area
        fx_first -= p[cell1] * normal[0] * area
        fx_linear -= p_linear * normal[0] * area

    # rho = 1, U = 1, A = D * span = 1  ->  C = 2 F  (eq. 9.98)
    return 2.0 * fx_first, 2.0 * fx_linear

if __name__ == "__main__":
    case = sys.argv[1]
    cd_first, cd_linear = drag_coefficients(case)
    print(f"{case}")
    print(f"  Cd pressione, prima cella      : {cd_first:.8f}")
    print(f"  Cd pressione, estrapolata      : {cd_linear:.8f}")
    print(f"  differenza                     : {cd_linear - cd_first:+.8f}")
    print("  (il primo valore deve coincidere con cd_pressure nel CSV)")