#set page(paper: "a4", margin: 2.5cm)
#set text(lang: "pl", size: 11pt)

#let przedmiot = "Algorytmy macierzowe"
#let prowadzacy = "[Prowadzący]"
#let grupa = "2"
#let autorzy = (
  (imie: "Wiktor Augustyn", album: "[Nr albumu 1]"),
  (imie: "Gabriel Wermiński", album: "[Nr albumu 2]"),
)

#align(center)[
  #text(size: 12pt)[#przedmiot]
  #v(0.6em)
  #text(size: 20pt, weight: "bold")[Sprawozdanie – Ćwiczenie 1]
  #v(0.3em)
  #text(size: 15pt)[Mnożenie macierzy]
]

#v(1em)
#line(length: 100%)
#v(0.5em)

#grid(
  columns: (1fr, 1fr),
  gutter: 1em,
  [
    *Autorzy:* \
    #for a in autorzy [
      #a.imie  \
    ]
  ],
  align(right)[
    *Grupa:* #grupa \ 
    *Data:* #datetime(year: 2026, month: 10, day: 15).display("[day].[month].[year]")
  ],
)

#v(0.5em)
#line(length: 100%)
#v(1.5em)

#import "@preview/lovelace:0.3.0": *

= Specyfikacja sprzętowa

Pomiary wykonano na laptopie Dell Inspiron 14 5420.

- Procesor: 12th Gen Intel(R) Core(TM) i5-1235U (4+8) @ 4.40 GHz.
- Pamięć operacyjna: 16 GB.
- System: EndeavourOS, Linux 7.2.9-arch1-1

= Opis algorytmów

W obu algorytmach zakładamy, że macierze są kwadratowe o rozmiarze $n times n$, gdzie $n$ jest potęgą dwójki (w przeciwnym razie macierze uzupełniamy zerami). Dla $n <= n_0$ stosujemy klasyczne mnożenie.

== Rekurencyjne mnożenie metodą Bineta

#figure(
  kind: "algorithm",
  supplement: [Algorytm],
  caption: [Rekurencyjne mnożenie macierzy metodą Bineta],
  pseudocode-list[
    + *funkcja* Binet($A$, $B$, $n$)
      + *jeśli* $n <= n_0$ *to zwróć* $A dot B$ (mnożenie klasyczne)
      + $k <- n slash 2$
      + podziel $A$ na bloki $A_(11), A_(12), A_(21), A_(22)$ o rozmiarze $k times k$
      + podziel $B$ na bloki $B_(11), B_(12), B_(21), B_(22)$ o rozmiarze $k times k$
      + $C_(11) <-$ Binet($A_(11)$, $B_(11)$, $k$) $+$ Binet($A_(12)$, $B_(21)$, $k$)
      + $C_(12) <-$ Binet($A_(11)$, $B_(12)$, $k$) $+$ Binet($A_(12)$, $B_(22)$, $k$)
      + $C_(21) <-$ Binet($A_(21)$, $B_(11)$, $k$) $+$ Binet($A_(22)$, $B_(21)$, $k$)
      + $C_(22) <-$ Binet($A_(21)$, $B_(12)$, $k$) $+$ Binet($A_(22)$, $B_(22)$, $k$)
      + złóż $C$ z bloków $C_(11), C_(12), C_(21), C_(22)$
      + *zwróć* $C$
  ],
)

Algorytm wykonuje 8 rekurencyjnych mnożeń macierzy rozmiaru $n slash 2$, więc $T(n) = 8 T(n slash 2) + Theta(n^2)$, czyli złożoność wynosi $Theta(n^3)$.

== Rekurencyjne mnożenie metodą Strassena

#figure(
  kind: "algorithm",
  supplement: [Algorytm],
  caption: [Rekurencyjne mnożenie macierzy metodą Strassena],
  pseudocode-list[
    + *funkcja* Strassen($A$, $B$, $n$)
      + *jeśli* $n <= n_0$ *to zwróć* $A dot B$ (mnożenie klasyczne)
      + $k <- n slash 2$
      + podziel $A$ na bloki $A_(11), A_(12), A_(21), A_(22)$ o rozmiarze $k times k$
      + podziel $B$ na bloki $B_(11), B_(12), B_(21), B_(22)$ o rozmiarze $k times k$
      + $M_1 <-$ Strassen($A_(11) + A_(22)$, $B_(11) + B_(22)$, $k$)
      + $M_2 <-$ Strassen($A_(21) + A_(22)$, $B_(11)$, $k$)
      + $M_3 <-$ Strassen($A_(11)$, $B_(12) - B_(22)$, $k$)
      + $M_4 <-$ Strassen($A_(22)$, $B_(21) - B_(11)$, $k$)
      + $M_5 <-$ Strassen($A_(11) + A_(12)$, $B_(22)$, $k$)
      + $M_6 <-$ Strassen($A_(21) - A_(11)$, $B_(11) + B_(12)$, $k$)
      + $M_7 <-$ Strassen($A_(12) - A_(22)$, $B_(21) + B_(22)$, $k$)
      + $C_(11) <- M_1 + M_4 - M_5 + M_7$
      + $C_(12) <- M_3 + M_5$
      + $C_(21) <- M_2 + M_4$
      + $C_(22) <- M_1 - M_2 + M_3 + M_6$
      + złóż $C$ z bloków $C_(11), C_(12), C_(21), C_(22)$
      + *zwróć* $C$
  ],
)

Algorytm wykonuje 7 rekurencyjnych mnożeń zamiast 8 (kosztem większej liczby dodawań), więc $T(n) = 7 T(n slash 2) + Theta(n^2)$, czyli złożoność wynosi $Theta(n^(log_2 7)) approx Theta(n^(2{,}81))$.

Na wykresach w następnym rozdziale linia przerywana i kropkowana to eksperymentalne dopasowanie prostej i paraboli do zmierzonych wyników Bineta i Strassena. To tylko trend. Oszacowaniem złożoności pozostają rekurencje powyżej, bo liczba operacji rośnie szybciej niż kwadrat.

= Fragmenty implementacji

Oba algorytmy dopełniają wejście zerami do najbliższej potęgi dwójki, a rekurencja schodzi do jednego mnożenia skalarnego.

Dla Bineta każdy krok rozcina macierz na cztery bloki i wykonuje osiem mnożeń rekurencyjnych. Dodanie dwóch bloków jest liczone jako tyle operacji zmiennoprzecinkowych, ile ma elementów wynik.

```python
n = next_pow2(max(m_A, n_A, m_B, n_B))
A_padded = np.pad(A, ((0, n - m_A), (0, n - n_A)), mode="constant")
B_padded = np.pad(B, ((0, n - m_B), (0, n - n_B)), mode="constant")

if A_i_e - A_i_s == 1:
    m_ops += 1
    return np.array([[A_padded[A_i_s][A_j_s] * B_padded[B_i_s][B_j_s]]])

C_11 = (rec(*se_11_A[0], *se_11_A[1], *se_11_B[0], *se_11_B[1]) +
        rec(*se_12_A[0], *se_12_A[1], *se_21_B[0], *se_21_B[1]))
a_s_ops += C_11.size
C_12 = (rec(*se_11_A[0], *se_11_A[1], *se_12_B[0], *se_12_B[1]) +
        rec(*se_12_A[0], *se_12_A[1], *se_22_B[0], *se_22_B[1]))
a_s_ops += C_12.size
C_21 = (rec(*se_21_A[0], *se_21_A[1], *se_11_B[0], *se_11_B[1]) +
        rec(*se_22_A[0], *se_22_A[1], *se_21_B[0], *se_21_B[1]))
a_s_ops += C_21.size
C_22 = (rec(*se_21_A[0], *se_21_A[1], *se_12_B[0], *se_12_B[1]) +
        rec(*se_22_A[0], *se_22_A[1], *se_22_B[0], *se_22_B[1]))
a_s_ops += C_22.size
return np.block([[C_11, C_12], [C_21, C_22]])
```

Strassen w każdym kroku liczy siedem iloczynów pomocniczych, a dodawania i odejmowania bloków też wchodzą do licznika operacji.

```python
M1 = rec(add(A11, A22), add(B11, B22))
M2 = rec(add(A21, A22), B11)
M3 = rec(A11, sub(B12, B22))
M4 = rec(A22, sub(B21, B11))
M5 = rec(add(A11, A12), B22)
M6 = rec(sub(A21, A11), add(B11, B12))
M7 = rec(sub(A12, A22), add(B21, B22))

C11 = add(sub(add(M1, M4), M5), M7)
C12 = add(M3, M5)
C21 = add(M2, M4)
C22 = add(sub(add(M1, M3), M2), M6)
```

// Ustaw na true, gdy pliki .png już istnieją w projekcie
#let pokaz-obrazy = true

#let wykres(plik, opis) = figure(
  if pokaz-obrazy {
    image(plik, width: 90%)
  } else {
    rect(width: 90%, height: 7cm, stroke: (dash: "dashed"), fill: luma(240))[
      #align(center + horizon)[`#plik`]
    ]
  },
  caption: opis,
)

= Wyniki pomiarów

Zmierzono rozmiary od 2 do 100. Rozmiar 1 jest przypadkiem bazowym: jedno mnożenie skalarne. Rozmiar 100 jest największym, który ta implementacja w czystym Pythonie kończy w akceptowalnym czasie. Dla $n = 1000$ macierz jest dopełniana do 1024, a sam Binet wykonuje wtedy $1024^3$ mnożeń skalarnych w interpreterze, więc ten rozmiar nie jest praktycznie obliczalny.

Operacje zmiennoprzecinkowe liczymy jako sumę mnożeń oraz dodawań i odejmowań. AlphaTensor występuje tylko dla rozmiarów, dla których plik rozkładów ma algorytm: 2, 3, 4, 5, 9, 10 i 11. Próbka jest za mała na dopasowanie krzywej.

== Czas wykonania w zależności od rozmiaru macierzy

#wykres(
  "imgs/multiplication_time.png",
  [Czas wykonania w zależności od rozmiaru macierzy $n$],
)

Przy $n = 100$ Strassen jest wyraźnie szybszy od Bineta. Skoki czasu wynikają z dopełniania do potęgi dwójki: na przykład 50 i 100 są liczone jako 64 i 128.

== Liczba operacji mnożenia w zależności od rozmiaru macierzy

#wykres(
  "imgs/multiplications_count.png",
  [Liczba operacji mnożenia w zależności od rozmiaru macierzy $n$],
)

Dla $n = 100$ (po dopełnieniu do 128) Binet wykonuje $128^3 = 2 thin 097 thin 152$ mnożeń, a Strassen $7^7 = 823 thin 543$.

== Liczba operacji dodawania i odejmowania w zależności od rozmiaru macierzy

#wykres(
  "imgs/numberOfAdditionsAndSubs.png",
  [Liczba operacji dodawania i odejmowania w zależności od rozmiaru macierzy $n$],
)

Krzywa Strassena leży wyraźnie wyżej. To koszt siedmiu iloczynów złożonych z sum i różnic bloków. Operacje zmiennoprzecinkowe to suma tych dwóch wykresów: Strassen ma ich łącznie więcej niż Binet, bo oszczędność na mnożeniach nie równoważy dodatkowych dodawań. Mimo to jest szybszy, bo mnożenia skalarne w tej implementacji ważą więcej niż dodawania.

= Kompilacja

Żeby program nie został w postaci interpretowanego skryptu Pythona, notatnik zamieniliśmy na plik źródłowy, a następnie skompilowaliśmy go narzędziem Nuitka. Nuitka tłumaczy ten plik na C i buduje z niego samodzielny plik wykonywalny. Do linkowania na Linuksie potrzebny jest `patchelf`.

```bash
jupyter nbconvert --to script lab1.ipynb --output lab1.py
nuitka --standalone lab1.py
```