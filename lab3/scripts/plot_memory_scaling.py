"""Wykres skalowania zuzycia pamieci w funkcji wymiaru n -> typst/imgs/memory_scaling_n.png

Uruchamiac z katalogu lab3:  python3 scripts/plot_memory_scaling.py
Dane: results/results.json (pole "scaling").
Pamiec = liczba zapisanych liczb f64 * 8 B * 3 kanaly (RGB), w MiB.
"""
import json

import matplotlib.pyplot as plt

BYTES_F64 = 8
CHANNELS = 3
MIB = 2**20

with open("results/results.json") as f:
    scaling = json.load(f)["scaling"]

ns = [s["n"] for s in scaling]
dense = [n * n * BYTES_F64 * CHANNELS / MIB for n in ns]
r1 = [s["hmat_r1_stored_floats"] * BYTES_F64 * CHANNELS / MIB for s in scaling]
r4 = [s["hmat_r4_stored_floats"] * BYTES_F64 * CHANNELS / MIB for s in scaling]

fig, ax = plt.subplots(figsize=(10, 5.6))
ax.plot(ns, dense, marker="s", color="#d62728", lw=2.2,
        label=r"Macierz gęsta ($3n^2$ liczb f64)")
ax.plot(ns, r1, marker="o", color="#1f77b4", lw=2.2,
        label=r"H-macierz ($r=1$, $\delta=\sigma_{n/2}$)")
ax.plot(ns, r4, marker="^", color="#2ca02c", lw=2.2,
        label=r"H-macierz ($r=4$, $\delta=3\,\sigma_{n/2}$)")

ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xticks(ns)
ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
ax.minorticks_off()
ax.set_xlabel(r"Rozmiar boku macierzy $n$ [piksele]", fontsize=12)
ax.set_ylabel("Zapotrzebowanie na pamięć [MiB] (skala log), RGB", fontsize=12)
ax.set_title(r"Skalowanie rozmiaru pamięci w funkcji wymiaru $n$",
             fontsize=14, fontweight="bold")
ax.grid(True, linestyle="--", alpha=0.5)
ax.legend(loc="upper left", fontsize=11)

fig.tight_layout()
fig.savefig("typst/imgs/memory_scaling_n.png", dpi=200)
print("zapisano typst/imgs/memory_scaling_n.png")
