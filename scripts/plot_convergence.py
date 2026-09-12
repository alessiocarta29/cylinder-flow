#!/usr/bin/env python3
# Grafico della convergenza di griglia: differenze fra livelli successivi
# in scala doppio logaritmica. La pendenza e' l'ordine osservato.
#
# Uso:  python3 scripts/plot_convergence.py

import matplotlib
matplotlib.use("Agg")          # salva su file senza aprire una finestra
import matplotlib.pyplot as plt

# Cd per livello, dai risultati dello studio di griglia
pressure_cell   = [1.17867501, 1.21578860, 1.24041912, 1.25265807, 1.25810840]
pressure_extrap = [1.21288056, 1.24340124, 1.25954282, 1.26347363, 1.26377081]
viscous         = [0.90131684, 0.87398654, 0.85224378, 0.84031633, 0.83408741]

# passo relativo alla griglia piu' grossa, associato alla prima delle due griglie
h = [1.0, 0.5, 0.25, 0.125]

def differences(c):
    return [abs(c[i + 1] - c[i]) for i in range(4)]

fig, ax = plt.subplots(figsize=(6.5, 5))
ax.loglog(h, differences(pressure_cell),   "o-", label="Pressure, cell value")
ax.loglog(h, differences(pressure_extrap), "s-", label="Pressure, extrapolated")
ax.loglog(h, differences(viscous),         "^-", label="Viscous")

# rette di riferimento, ancorate ciascuna alla curva che descrive
a1 = differences(viscous)[-1]
ax.loglog(h, [a1 * x / h[-1] for x in h], "k--", lw=0.8, label="first order")
a2 = differences(pressure_extrap)[2]
ax.loglog(h, [a2 * (x / h[2]) ** 2 for x in h], "k:", lw=0.8, label="second order")

ax.set_xticks(h)
ax.set_xticklabels(["1", "1/2", "1/4", "1/8"])
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())   # via le etichette intermedie
ax.set_xlabel("h relative to the coarsest grid")
ax.set_ylabel("|difference| in Cd between successive levels")
ax.set_title("Grid convergence, cylinder at Re = 20")
ax.grid(True, which="both", ls=":", lw=0.5)
ax.legend()
fig.tight_layout()
fig.savefig("figures/grid_convergence_re20.png", dpi=200)
print("scritto figures/grid_convergence_re20.png")