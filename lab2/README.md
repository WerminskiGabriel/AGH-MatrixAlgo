# Lab 2 — odwracanie macierzy, faktoryzacja LU i algorytm Ladermana

Rekurencyjne odwracanie macierzy oraz faktoryzacja LU dla algorytmów Bineta, Strassena i Ladermana, a także rekurencyjne mnożenie macierzy metodą Ladermana. Są pomiary czasu i liczby operacji (dla $n = 1..100$ oraz $n = 1..1000$) oraz sprawozdanie w Typst.

Komendy poniżej odpalaj z katalogu `lab2`.

## Kompilacja Cython

Algorytmy w `matrix_algos_cy.pyx` są kompilowane przez Cythona i GCC do natywnego modułu C (`.so`). Notatnik sam uruchamia kompilację przy pierwszym starcie, ale można też zbudować moduł ręcznie:

```bash
python setup.py build_ext --inplace
```

## Notatnik

```bash
jupyter lab lab2.ipynb
```

## Sprawozdanie

```bash
typst compile typst/main.typ typst/output.pdf
```
