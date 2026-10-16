# Lab 1 — mnożenie macierzy

Porównanie trzech algorytmów mnożenia macierzy kwadratowych: Bineta, Strassena i AlphaTensor. Są pomiary czasu i liczby operacji oraz sprawozdanie w Typst.

Komendy poniżej odpalaj z katalogu `lab1`.

## Notatnik

```bash
jupyter lab lab1.ipynb
```

## Sprawozdanie

```bash
typst compile typst/main.typ typst/output.pdf
```

## Kompilacja

`patchelf` jest potrzebny, żeby Nuitka zlinkowała program.

```bash
jupyter nbconvert --to script lab1.ipynb --output lab1.py
nuitka --standalone lab1.py
./lab1.dist/lab1.bin
```
