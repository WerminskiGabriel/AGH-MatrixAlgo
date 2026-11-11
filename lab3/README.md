# Lab 3 — Hierarchiczna kompresja macierzy (H-macierze) i zredukowane SVD (rSVD)

Projekt z przedmiotu *Algorytmy macierzowe* (AGH, semestr 5).
**Autorzy:** Wiktor Augustyn, Gabriel Wermiński (Grupa 2).

Implementacja rekurencyjnej kompresji macierzy reprezentujących 3 kanały barwne ($R, G, B \in [0, 255]^{512 \times 512}$) kolorowej bitmapy wysokiej rozdzielczości ($512 \times 512$, $2^9 \times 2^9$) w oparciu o hierarchiczne macierze (H-macierze, quadtree) oraz zredukowaną dekompozycję według wartości osobliwych (*truncated/thin SVD*).

---

## Architektura rozwiązania

1. **Natywny moduł numeryczny w języku Rust (`src/lib.rs`, `src/main.rs`, `Cargo.toml`):**
   - Skompilowany z flagami maksymalnej optymalizacji (`opt-level = 3`, `lto = "fat"`, `codegen-units = 1`) oraz wykorzystaniem instrukcji wektorowych AVX2/FMA za pośrednictwem biblioteki `faer`.
   - Zapewnia brak narzutu interpretera i Garbage Collectora, deterministyczny układ pamięci drzewa czwórkowego (`Box<[HNode; 4]>`) oraz wydajność rzędu $20 - 53\text{ GFLOP/s}$ (do $150\%$ mocy szczytowej pojedynczego rdzenia CPU).
   - Realizuje algorytmy: `create_tree`, `compress_matrix`, `decompress_hmatrix`, szybkie mnożenie H-macierzy przez wektor (`hmatrix_matvec`), zliczanie operacji zmiennoprzecinkowych (FLOPs), rysowanie struktury blokowej H-macierzy z pasami $U, V$ oraz eksport bitmap i metryk numerycznych do `results/results.json` i `results/small_test.json`.
2. **Notatnik Jupyter (`lab3.ipynb`):**
   - Kernel Python 3, automatyczna orkiestracja kompilacji i uruchomienia binarnego modułu Rust, weryfikacja numeryczna z pakietem **NumPy** oraz **GNU Octave**, a także generowanie wykresów naukowej jakości (`matplotlib`) do katalogu `typst/imgs/`.
3. **Weryfikacja referencyjna w GNU Octave (`scripts/verify_octave.m`):**
   - Niezależny skrypt GNU Octave / MATLAB weryfikujący zgodność wartości osobliwych i błędu rekonstrukcji do precyzji maszynowej ($< 10^{-14}$).
4. **Sprawozdanie w Typst (`typst/main.typ` $\to$ `typst/output.pdf`):**
   - Pełne sprawozdanie (18 stron) zachowujące układ graficzny, szatę typograficzną, makra i tabele zgodne z `lab2`.

---

## Struktura plików w katalogu `lab3/`

```text
lab3/
├── Cargo.toml                  # Konfiguracja projektu Rust (faer, image, serde, serde_json, rayon)
├── src/
│   ├── lib.rs                  # Główny moduł biblioteki i re-eksporty API
│   ├── tree.rs                 # Struktura HNode, algorytmy create_tree, compress_matrix, decompress, matvec
│   ├── svd_utils.rs            # Wrappery thin-SVD, estymacja FLOPs, metryki błędu Frobeniusa i PSNR
│   ├── visualizer.rs           # Operacje na bitmapach, zapis kanałów RGB oraz rysowacz struktury H-macierzy
│   ├── io_stats.rs             # Struktury wyników eksperymentów oraz eksport do formatów JSON i CSV
│   └── main.rs                 # Główna funkcja sterująca przebiegiem eksperymentów (CLI / potok pomiarowy)
├── scripts/
│   ├── verify_octave.m         # Skrypt GNU Octave do weryfikacji numerycznej SVD
├── assets/
│   └── input_512.png           # Wejściowa kolorowa bitmapa 512x512 (RGB)
├── results/
│   ├── results.json            # Zrzut wszystkich pomiarów, wartości osobliwych i statystyk dla Pythona
│   └── small_test.json         # Dane małych macierzy (4x4, 16x16) do testów numerycznych z NumPy i Octave
├── typst/
│   ├── imgs/                   # Wygenerowane wykresy matplotlib i bitmapy kompresji (.png)
│   ├── main.typ                # Kod źródłowy sprawozdania w Typst
│   └── output.pdf              # Skompilowane sprawozdanie PDF
├── lab3.ipynb                  # Wykonany notatnik Jupyter z zapisanymi outputami i wykresami
└── README.md                   # Niniejsza instrukcja
```

---

## Instrukcja uruchomienia

Wszystkie poniższe polecenia należy wykonywać będąc w katalogu `lab3/`.

### 1. Kompilacja i uruchomienie modułu w Rust

Uruchomienie profilu wydaniowego (`--release`):

```bash
cargo run --release
```

Program automatycznie:
- Wczyta bitmapę wejściową `assets/input_512.png` i rozdzieli ją na macierze kanałów $R, G, B \in [0, 255]^{512 \times 512}$.
- Wyznaczy pełny rozkład SVD dla każdego kanału i znajdzie wartości osobliwe $\sigma_1, \sigma_{256}, \sigma_{512}$.
- Przeprowadzi kompresję dla 6 zadanych konfiguracji oraz 2 konfiguracji optymalnych (`case1`..`case6`, `opt1`, `opt2`).
- Wygeneruje rysunki struktur H-macierzy (`{case}_tree_{r,g,b}.png`), zrekonstruowane kanały oraz obrazy RGB (`{case}_rec_rgb.png`).
- Przeprowadzi sweep parametrów $(r, \delta)$ dla krzywych Rate-Distortion oraz skalowanie $N \in [32, 1024]$ dla wydajności GFLOP/s i `% Peak FLOP/s`.
- Zapisze wyniki do `results/results.json` oraz `results/small_test.json`.

### 2. Notatnik Jupyter i wykresy

Uruchomienie i przeliczenie notatnika (w środowisku wirtualnym `.venv`):

```bash
# Interaktywne otwarcie w przeglądarce:
jupyter lab lab3.ipynb

# Lub automatyczne wykonanie wsadowe w miejscu:
jupyter nbconvert --to notebook --execute --inplace lab3.ipynb
```

### 3. Weryfikacja w GNU Octave

Uruchomienie skryptu weryfikacyjnego GNU Octave:

```bash
octave --no-gui --quiet scripts/verify_octave.m
```

### 4. Kompilacja sprawozdania Typst

Kompilacja pliku `main.typ` do `output.pdf`:

```bash
typst compile typst/main.typ typst/output.pdf
```

---

## Podsumowanie wyników

| Konfiguracja | Rząd $r$ | Próg $\delta$ | Stopień kompresji (CR) | RAM [KB] | Oszczędność RAM | Jakość PSNR [dB] | Czas [ms] | Wydajność (% Peak CPU) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Case 1** | 1 | $\sigma_1$ | **255.8x** | 24.4 | 99.6% | 17.25 | 159.2 | 53.11 GFLOP/s (150.9%) |
| **Case 2** | 1 | $\sigma_{512}$ | **0.81x** | 64 216.0 | -23.7% | 39.22 | 1195.1 | 14.12 GFLOP/s (40.1%) |
| **Case 3** | 1 | $\sigma_{256}$ | **3.49x** | 8 645.5 | 71.3% | 35.40 | 913.8 | 17.94 GFLOP/s (51.0%) |
| **Case 4** | 4 | $\sigma_1$ | **63.9x** | 96.4 | 98.4% | 20.68 | 172.6 | 49.00 GFLOP/s (139.2%) |
| **Case 5** | 4 | $\sigma_{512}$ | **0.45x** | 27 858.6 | -122.9% | 99.00 | 1154.2 | 14.59 GFLOP/s (41.5%) |
| **Case 6** | 4 | $\sigma_{256}$ | **2.30x** | 3 674.6 | 56.5% | 43.30 | 776.0 | 20.79 GFLOP/s (59.1%) |
| **Opt 1** | 2 | $3\sigma_{256}$ | **4.89x** | 2 105.8 | **79.6%** | **31.76** | 789.4 | 20.43 GFLOP/s (58.0%) |
| **Opt 2** | 4 | $4\sigma_{256}$ | **4.39x** | 1 640.8 | **77.2%** | **32.82** | 706.6 | 22.60 GFLOP/s (64.2%) |
