"""Wykres czasu i liczby FLOPs w funkcji wymiaru n -> typst/imgs/time_and_flops_n.png

Uruchamiac z katalogu lab3:  python3 scripts/plot_time_flops.py
Dane: results/results.json (pole "scaling").
"""
import json

import matplotlib.pyplot as plt

with open("results/results.json") as f:
    scaling = json.load(f)["scaling"]

ns = [s["n"] for s in scaling]
series = [
    ("Pełny rozkład SVD", "full_svd", "#d62728", "s"),
    (r"H-macierz ($r=1$)", "hmat_r1", "#1f77b4", "o"),
    (r"H-macierz ($r=4$)", "hmat_r4", "#2ca02c", "^"),
]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.0))

for label, key, color, marker in series:
    ax1.plot(ns, [s[f"{key}_time_s"] for s in scaling], marker=marker, color=color, lw=2.0, label=label)
    ax2.plot(ns, [s[f"{key}_flops"] for s in scaling], marker=marker, color=color, lw=2.0, label=label)

ax1.set_title(r"Czas wykonania w funkcji wymiaru $n$", fontsize=12, fontweight="bold")
ax1.set_ylabel("Czas wykonania [s] (skala log)", fontsize=11)
ax2.set_title("Liczba operacji zmiennoprzecinkowych (FLOPs)", fontsize=12, fontweight="bold")
ax2.set_ylabel("Liczba operacji zmiennoprzecinkowych [FLOPs]", fontsize=11)

for ax in (ax1, ax2):
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(ns)
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.minorticks_off()
    ax.set_xlabel(r"Wymiar macierzy $n$", fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="upper left", fontsize=10)

fig.tight_layout()
fig.savefig("typst/imgs/time_and_flops_n.png", dpi=200)
print("zapisano typst/imgs/time_and_flops_n.png")
