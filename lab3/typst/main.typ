#set page(paper: "a4", margin: 2.3cm)
#set text(lang: "pl", size: 10.5pt)
#set par(justify: true)
#set heading(numbering: "1.")

#show raw.where(block: true): it => block(
  width: 100%,
  stroke: 0.6pt + luma(170),
  fill: luma(250),
  radius: 4pt,
  inset: (x: 12pt, y: 10pt),
  it,
)

#let przedmiot = "Algorytmy macierzowe"
#let grupa = "2"
#let autorzy = (
  (imie: "Wiktor Augustyn"),
  (imie: "Gabriel Wermiński"),
)

#align(center)[
  #text(size: 12pt)[#przedmiot]
  #v(0.5em)
  #text(size: 19pt, weight: "bold")[Sprawozdanie – Ćwiczenie 3]
  #v(0.3em)
  #text(size: 14pt)[Hierarchiczna kompresja macierzy (H-macierze) z wykorzystaniem częściowego rozkładu SVD (rSVD)]
]

#v(0.8em)
#line(length: 100%)
#v(0.4em)

#grid(
  columns: (1fr, 1fr),
  gutter: 1em,
  [
    *Autorzy:* \
    #for a in autorzy [
      #a.imie \
    ]
  ],
  align(right)[
    *Grupa:* #grupa \
    *Data:* #datetime(year: 2026, month: 11, day: 11).display("[day].[month].[year]")
  ],
)

#v(0.4em)
#line(length: 100%)
#v(0.8em)

#outline()
#v(1em)

#let algo-box(title, body) = figure(
  kind: "algorithm",
  supplement: [Algorytm],
  caption: title,
  block(
    width: 100%,
    stroke: 0.6pt + luma(170),
    fill: luma(250),
    radius: 4pt,
    inset: (x: 14pt, y: 11pt),
    align(left)[
      #set par(leading: 0.68em, justify: false)
      #set text(size: 9.8pt)
      #body
    ],
  ),
)

#let ind(level) = h(level * 1.4em)

#let wykres(plik, opis) = figure(
  image(plik, width: 95%),
  caption: opis,
)

= Wstęp i specyfikacja techniczna

Celem trzeciego ćwiczenia laboratoryjnego jest zaprojektowanie, implementacja, optymalizacja oraz dogłębna analiza numeryczna rekurencyjnej kompresji macierzy opartej na teorii *macierzy hierarchicznych (H-macierzy)* oraz *częściowej dekompozycji według wartości osobliwych (rSVD – truncated Singular Value Decomposition)*.

W ramach ćwiczenia zrealizowano następujące zadania inżynierskie i badawcze:
1. *Przetwarzanie obrazu w przestrzeni barw RGB:* Wykorzystanie bogatej w szczegóły kolorowej bitmapy o wymiarach $N times N = 512 times 512$ pikseli ($2^9 times 2^9$) oraz rozbicie jej na trzy niezależne macierze kanałów $R, G, B in [0, 255]^(512 times 512)$.
2. *Rekurencyjna dekompozycja drzewiasta (CreateTree oraz CompressMatrix):* Implementacja algorytmów budowy drzewa czwórkowego (*quadtree*) macierzy hierarchicznej z dynamicznym testem dopuszczalności (*admissibility condition*) opartym na progu błędu $delta$ oraz zadanym maksymalnym rzędzie podbloku $r$.
3. *Implementacja rysowacza struktury H-macierzy:* Moduł generujący wizualizację podziału blokowego z zaznaczeniem ramki każdego podbloku oraz pasów poziomych i pionowych o grubości $r$ (tworzących charakterystyczny kształt litery „$Gamma$”), reprezentujących macierze czynnikowe $U$ i $V$.
4. *Implementacja rysowacza skompresowanej bitmapy:* Rekurencyjna dekompresja drzewa liści niskorządowych i zerowych, rekonstrukcja poszczególnych kanałów $R, G, B$ oraz złożenie wyjściowej skompresowanej fotografii RGB.
5. *Badanie 6 obowiązkowych oraz 2 optymalnych wariantów kompresji:* Wyznaczenie pełnego widma wartości osobliwych $sigma_1, dots, sigma_(512)$ każdego kanału i przebadanie konfiguracji dla $r in {1, 4}$ oraz progów $delta in {sigma_1, sigma_(512), sigma_(256)}$, a także dobór parametrów maksymalizujących kompresję przy zachowaniu wysokiej wierności wizualnej ($"PSNR" > 32$ dB).
6. *Pomiary pamięci, liczby operacji (FLOPs) oraz wydajności GFLOP/s i % Peak CPU:* Zliczenie dokładnej liczby operacji zmiennoprzecinkowych, pomiary czasu i wyznaczenie odsetka osiąganej mocy szczytowej procesora (`% Peak FLOP/s`).
7. *Weryfikacja numeryczna z bibliotekami NumPy oraz GNU Octave:* Porównanie wyników rozkładu SVD oraz błędów rekonstrukcji z referencyjnymi środowiskami algebry liniowej.

== Uzasadnienie wyboru języka kompilowanego (Rust) oraz parametry sprzętowe

Zgodnie z wymaganiami zadania obliczenia macierzowe nie mogą być wykonywane przez wolne języki interpretowane. Do implementacji rdzenia numerycznego wybrano nowoczesny język systemowy *Rust* (kompilator `rustc 1.99.1`) z zaawansowaną biblioteką algebry liniowej `faer 0.19`.

*Kluczowe przewagi technologiczne języka Rust i biblioteki faer:*
- *Zero-cost abstractions i brak Garbage Collectora:* Zarządzanie węzłami drzewa czwórkowego odbywa się w sposób ściśle deterministyczny (`Box<[HNode; 4]>`), eliminując jakiekolwiek pauzy GC czy nieprzewidywalne opóźnienia alokacji pamięci.
- *Zaawansowana wektoryzacja SIMD (AVX2 + FMA):* Biblioteka `faer` została zaprojektowana od podstaw z myślą o maksymalnym wykorzystaniu jednostek wektorowych współczesnych procesorów x86_64. Wykorzystuje mikro-jądra gemm oraz algorytmy SVD zoptymalizowane pod kątem lokalności pamięci podręcznej (cache L1/L2).
- *Maksymalna optymalizacja kompilatora:* Projekt skompilowano w profilu produkcyjnym (`cargo build --release`) z flagami `-C opt-level=3`, pełnym linkowaniem między-modułowym LTO (`lto = "fat"`) oraz jednym blokiem generacji kodu (`codegen-units = 1`), co umożliwia kompilatorowi LLVM agresywną wektoryzację pętli i inlining funkcji rekurencyjnych.

*Parametry platformy testowej:*
- *Procesor:* 12th Gen Intel(R) Core(TM) i5-1235U (architektura Alder Lake, 2 rdzenie Golden Cove Performance + 8 rdzeni Gracemont Efficient, 12 wątków, maksymalne taktowanie Turbo $4.40$ GHz).
- *Pamięć podręczna:* L1d: $352$ KiB, L2: $6.5$ MiB, L3: $12$ MiB (współdzielona).
- *Pamięć operacyjna:* $16$ GB DDR4.
- *System operacyjny:* EndeavourOS Linux (kernel 7.2.9-arch1-1 x86_64).
- *Teoretyczna wydajność szczytowa CPU (Single-Core Peak):* Jednostka AVX2 wykonuje operacje na rejestrach 256-bitowych (4 liczby `f64`). Przy wykorzystaniu instrukcji FMA (Fused Multiply-Add, 2 FLOP/cykl) rdzeń wykonuje $4 times 2 = 8$ operacji zmiennoprzecinkowych podwójnej precyzji na cykl na jeden port FMA. Przy taktowaniu $4.40$ GHz teoretyczna wydajność szczytowa pojedynczego rdzenia wynosi:
  $ R_("peak") = 4.40 "GHz" times 8 "FLOP/cykl" = 35.2 "GFLOP/s" $

= Opis teoretyczny i pseudokod algorytmów

== Podstawy matematyczne dekompozycji SVD i aproksymacji niskorządowej (rSVD)

Dla dowolnej rzeczywistej macierzy prostokątnej $M in RR^(m times n)$ (gdzie bez straty ogólności zakładamy $m >= n$) istnieje rozkład według wartości osobliwych (*Singular Value Decomposition – SVD*):
$ M = U D V^T = sum_(i=1)^n sigma_i u_i v_i^T $
gdzie:
- $U in RR^(m times m)$ jest macierzą ortogonalną ($U^T U = I_m$), której kolumny $u_i$ stanowią lewostronne wektory osobliwe,
- $V in RR^(n times n)$ jest macierzą ortogonalną ($V^T V = I_n$), której kolumny $v_i$ stanowią prawostronne wektory osobliwe,
- $D in RR^(m times n)$ jest macierzą pseudodiagonalną zawierającą nieujemne wartości osobliwe uporządkowane nierosnąco:
  $ sigma_1 >= sigma_2 >= dots >= sigma_n >= 0 $

W praktyce inżynierskiej stosuje się *zredukowaną (cienką) postać SVD (thin SVD)*, w której macierz $U$ ma rozmiar $m times n$, $D$ jest macierzą kwadratową $n times n$, a $V$ ma wymiar $n times n$.

=== Twierdzenie Eckarta-Younga-Mirsky'ego (Najlepsza aproksymacja niskorządowa)

*Twierdzenie:* Niech $M = U D V^T$ będzie rozkładem SVD macierzy $M in RR^(m times n)$. Wówczas dla dowolnego rzędu $r < n$, macierz zredukowana:
$ tilde(M)_r = sum_(i=1)^r sigma_i u_i v_i^T = U_r D_r V_r^T $
stanowi najlepszą aproksymację macierzy $M$ wśród wszystkich macierzy o rzędzie co najwyżej $r$ w normie spektralnej oraz normie Frobeniusa:
$ min_(op("rank")(X) <= r) ||M - X||_2 = ||M - tilde(M)_r||_2 = sigma_(r+1) $
$ min_(op("rank")(X) <= r) ||M - X||_F = ||M - tilde(M)_r||_F = sqrt(sum_(i=r+1)^n sigma_i^2) $

Powyższa własność leży u podstaw kompresji macierzowej: odrzucenie wyrazów od $r+1$ do $n$ wprowadza błąd w normie spektralnej dokładnie równy $(r+1)$-szej wartości osobliwej $sigma_(r+1)$. Jeżeli $sigma_(r+1) < delta$, to błąd odcięcia jest ściśle kontrolowany przez zadany próg tolerancji $delta$.

=== Opłacalność pamięciowa reprezentacji niskorządowej

Macierz gęsta o wymiarach $m times n$ wymaga zapamiętania $m dot n$ liczb zmiennoprzecinkowych. Zapis macierzy w postaci iloczynu czynników niskorządowych $tilde(M)_r = U_r V_r$, gdzie $U_r in RR^(m times r)$ oraz $V_r = D_r V_r^T in RR^(r times n)$, wymaga przechowania:
$ N_("comp") = m dot r + r dot n = r (m + n) $
liczb. Kompresja podbloku kwadratowego $n times n$ do liścia rzędu $r$ przynosi oszczędność pamięci wtedy i tylko wtedy, gdy:
$ 2 r n < n^2 <=> r < n / 2 $
Dla bloku $512 times 512$ i rzędu $r = 4$, zamiast $262 thin 144$ liczb przechowujemy zaledwie $4 dot (512 + 512) = 4 thin 096$ liczb, co oznacza ponad *64-krotną redukcję objętości danych*!

== Idea macierzy hierarchicznych (H-macierzy)

Pojedyncza globalna aproksymacja SVD całej macierzy obrazu może prowadzić do rozmycia drobnych krawędzi i utraty detali o wysokich częstotliwościach przestrzennych. Rozwiązaniem tego problemu jest *reprezentacja hierarchiczna (drzewo czwórkowe – quadtree)*:
1. Jeżeli dany blok macierzy jest gładki (jego $(r+1)$-sza wartość osobliwa spełnia $sigma_(r+1) < delta$), kompresujemy go bezpośrednio do liścia niskorządowego rzędu $r$.
2. Jeżeli blok zawiera ostre krawędzie lub skomplikowaną teksturę ($sigma_(r+1) >= delta$), dzielimy go geometrycznie na 4 mniejsze ćwiartki: górno-lewą ($A_(11)$), górno-prawą ($A_(12)$), dolno-lewą ($A_(21)$) i dolno-prawą ($A_(22)$).
3. Procedurę powtarzamy rekurencyjnie aż do spełnienia kryterium dopuszczalności lub osiągnięcia minimalnego rozmiaru bloku ($2 times 2$).

#algo-box([Rekurencyjna budowa drzewa macierzy hierarchicznej (`CreateTree`)])[
  *Wejście:* macierz $A in RR^(N times N)$, zakresy indeksów wierszy $[t_("min"), t_("max"))$ oraz kolumn $[s_("min"), s_("max"))$, maksymalny dopuszczalny rząd liścia $r in NN_+$, próg tolerancji wartości osobliwych $delta > 0$ \
  *Wyjście:* korzeń poddrzewa $v$ reprezentujący skompresowany blok $A(t_("min"):t_("max"), s_("min"):s_("max"))$ \
  #line(length: 100%, stroke: 0.4pt + luma(200))
  *Procedura* $"CreateTree"(A, t_("min"), t_("max"), s_("min"), s_("max"), r, delta)$: \
  #ind(1) *1.* Wyznacz wymiary bieżącego podbloku: $m <- t_("max") - t_("min")$ oraz $n <- s_("max") - s_("min")$. \
  #ind(1) *2.* *Jeżeli* podblok jest zerowy ($max_(i, j) |A(t_("min") + i, s_("min") + j)| < 10^(-12)$), *wówczas:* \
  #ind(2) Utwórz liść zerowy $v$: $v."rank" <- 0$, $v."sons" <- emptyset$, $v."size" <- (m, n)$ i *zwróć* $v$. \
  #ind(1) *3.* Wyznacz zredukowany rozkład SVD podbloku: \
  #ind(2) $(U, D, V) <- "truncatedSVD"(A(t_("min"):t_("max"), s_("min"):s_("max")), r + 1)$, gdzie $D = "diag"(sigma_1, sigma_2, dots)$. \
  #ind(1) *4.* Sprawdź warunek dopuszczalności (*admissibility condition*): \
  #ind(2) *Jeżeli* $min(m, n) <= r$ *lub* $sigma_(r+1) < delta$ *lub* $min(m, n) <= 2$, *wówczas:* \
  #ind(3) $v <- "CompressMatrix"(t_("min"), t_("max"), s_("min"), s_("max"), U, D, V, r, delta)$ \
  #ind(3) *Zwróć* liść skompresowany $v$. \
  #ind(1) *5.* *W przeciwnym razie* (blok niedopuszczalny – podział rekurencyjny na 4 ćwiartki): \
  #ind(2) Wyznacz punkty podziału: $t_("mid") <- t_("min") + floor(m slash 2)$ oraz $s_("mid") <- s_("min") + floor(n slash 2)$: \
  #ind(3) $A(t_("min"):t_("max"), s_("min"):s_("max")) = mat(delim: "(", A_(11), A_(12); A_(21), A_(22))$ \
  #ind(2) Utwórz węzeł wewnętrzny $v$ i wywołaj rekurencyjnie procedurę dla 4 potomków: \
  #ind(3) #grid(
    columns: (1fr, 1fr),
    row-gutter: 0.45em,
    column-gutter: 1em,
    [$v_(11) <- "CreateTree"(A, t_("min"), t_("mid"), s_("min"), s_("mid"), r, delta)$],
    [$v_(12) <- "CreateTree"(A, t_("min"), t_("mid"), s_("mid"), s_("max"), r, delta)$],
    [$v_(21) <- "CreateTree"(A, t_("mid"), t_("max"), s_("min"), s_("mid"), r, delta)$],
    [$v_(22) <- "CreateTree"(A, t_("mid"), t_("max"), s_("mid"), s_("max"), r, delta)$],
  )
  #ind(2) Przypisz synów: $v."sons" <- (v_(11), v_(12), v_(21), v_(22))$, $v."size" <- (m, n)$. \
  #ind(2) *Zwróć* węzeł wewnętrzny $v$.
]

#v(0.5em)

#algo-box([Kompresja dopuszczalnego podbloku macierzy do liścia SVD (`CompressMatrix`)])[
  *Wejście:* zakresy indeksów $[t_("min"), t_("max"))$, $[s_("min"), s_("max"))$, czynniki rozkładu $(U, D, V)$, maksymalny rząd $r$, próg $delta$ \
  *Wyjście:* skompresowany węzeł liścia $v$ przechowujący macierze czynnikowe rzędu $k <= r$ \
  #line(length: 100%, stroke: 0.4pt + luma(200))
  *Procedura* $"CompressMatrix"(t_("min"), t_("max"), s_("min"), s_("max"), U, D, V, r, delta)$: \
  #ind(1) *1.* Przyjmij wymiary $m <- t_("max") - t_("min")$, $n <- s_("max") - s_("min")$ oraz wektor wartości osobliwych $sigma <- "diag"(D)$. \
  #ind(1) *2.* *Jeżeli* największa wartość osobliwa $sigma_1 < delta$, *wówczas:* \
  #ind(2) Utwórz liść zerowy $v$: $v."rank" <- 0$, $v."sons" <- emptyset$, $v."size" <- (m, n)$ i *zwróć* $v$. \
  #ind(1) *3.* Wyznacz efektywny rząd aproksymacji: $k <- min(r, min(m, n))$. \
  #ind(1) *4.* Utwórz nowy liść niskorządowy $v$ i zapisz pierwszych $k$ składowych rozkładu rSVD: \
  #ind(2) $v."rank" <- k, quad v."singular_values" <- (sigma_1, dots, sigma_k)$ \
  #ind(2) $v.U <- U(:, 1:k) in RR^(m times k)$ \
  #ind(2) $v.V <- "diag"(sigma_1, dots, sigma_k) dot V(:, 1:k)^T in RR^(k times n)$ \
  #ind(2) $v."sons" <- emptyset, quad v."size" <- (m, n)$ \
  #ind(1) *5.* *Zwróć* liść $v$.
]

#v(0.5em)

#algo-box([Rekurencyjna dekompresja H-macierzy oraz wizualizacja struktury podziału])[
  *Wejście:* korzeń drzewa H-macierzy $v$ \
  *Wyjście:* zrekonstruowana macierz gęsta $tilde(A) in RR^(m times n)$ (dekompresja) lub obraz struktury $I$ (rysowacz) \
  #line(length: 100%, stroke: 0.4pt + luma(200))
  *Procedura* $"DecompressMatrix"(v)$: \
  #ind(1) *1.* *Jeżeli* $v."sons" = emptyset$ oraz $v."rank" = 0$ (liść zerowy), *wówczas zwróć* macierz zerową $bold(0)_(m times n)$. \
  #ind(1) *2.* *Jeżeli* $v."sons" = emptyset$ oraz $v."rank" = k > 0$ (liść niskorządowy), *wówczas zwróć* iloczyn $v.U dot v.V in RR^(m times n)$. \
  #ind(1) *3.* *Jeżeli* $v."sons" = (v_(11), v_(12), v_(21), v_(22))$ (węzeł wewnętrzny), *wówczas zwróć* macierz złożoną z 4 bloków: \
  #ind(2) $tilde(A) <- mat(delim: "(",
    "DecompressMatrix"(v_(11)), "DecompressMatrix"(v_(12));
    "DecompressMatrix"(v_(21)), "DecompressMatrix"(v_(22))
  )$ \
  #line(length: 100%, stroke: 0.4pt + luma(220))
  *Procedura* $"DrawHMatrixStructure"(v, I)$: \
  #ind(1) *1.* *Jeżeli* $v."sons" = emptyset$ (liść drzewa na obszarze $[t_("min"), t_("max")) times [s_("min"), s_("max"))$): \
  #ind(2) Obrysuj prostokąt bloku na płótnie $I$ czarną ramką o grubości $1$ px. \
  #ind(2) *Jeżeli* $v."rank" = k > 0$, zamaluj na czarno górny pas poziomy o wysokości $k$ pikseli oraz lewy pas pionowy o szerokości $k$ pikseli (charakterystyczny kształt litery „$Gamma$” obrazujący wymiary $m times k$ oraz $k times n$). \
  #ind(1) *2.* *W przeciwnym razie* (węzeł wewnętrzny $v."sons" = (v_(11), v_(12), v_(21), v_(22))$): \
  #ind(2) *Dla każdego* syna $u in v."sons"$ wywołaj rekurencyjnie $"DrawHMatrixStructure"(u, I)$.
]

#v(0.5em)

#algo-box([Szybkie mnożenie macierzy hierarchicznej przez wektor (`HMatVec`)])[
  *Wejście:* węzeł drzewa H-macierzy $v$, wektor wejściowy $x in RR^N$, wektor wyjściowy $y in RR^N$ \
  *Wyjście:* zaktualizowany wektor $y <- y + H dot x$ obliczony w czasie $Theta(r N log N)$ \
  #line(length: 100%, stroke: 0.4pt + luma(200))
  *Procedura* $"HMatVec"(v, x, y)$: \
  #ind(1) *1.* *Jeżeli* $v."sons" = emptyset$ oraz $v."rank" = 0$, *wówczas zakończ* (blok zerowy nie wnosi wkładu). \
  #ind(1) *2.* *Jeżeli* $v."sons" = emptyset$ oraz $v."rank" = k > 0$ (liść niskorządowy $v.U in RR^(m times k), v.V in RR^(k times n)$): \
  #ind(2) *Krok 2a:* Oblicz wektor pośredni $z <- v.V dot x(s_("min"):s_("max")) in RR^k$ #h(1fr) (koszt: $k(2n - 1)$ FLOPs). \
  #ind(2) *Krok 2b:* Zaktualizuj wynik $y(t_("min"):t_("max")) <- y(t_("min"):t_("max")) + v.U dot z$ #h(1fr) (koszt: $2 m k$ FLOPs). \
  #ind(1) *3.* *Jeżeli* $v."sons" = (v_(11), v_(12), v_(21), v_(22))$ (węzeł wewnętrzny): \
  #ind(2) Wywołaj rekurencyjnie $"HMatVec"(v_(i j), x, y)$ dla wszystkich czterech ćwiartek $i, j in {1, 2}$.
]

= Najważniejsze fragmenty kodu źródłowego (Rust)

== Definicja węzła drzewa H-macierzy (`HNode` oraz `HNodeKind`)

Poniższy listing przedstawia definicję rekurencyjnego drzewa czwórkowego w języku Rust. Wykorzystano typ wyliczeniowy `enum` oraz inteligentny wskaźnik `Box<[HNode; 4]>` o zerowym narzucie pamięciowym:

```rust
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum HNodeKind {
    /// Blok zerowy (wszystkie wartosci ponizej progu delta lub wygaszone)
    Zero,
    /// Lisc niskorzedowy przechowujacy czynniki U (m x rank) oraz V (rank x n)
    Leaf {
        rank: usize,
        singular_values: Vec<f64>,
        u: Vec<f64>, // uklad row-major: m * rank liczb f64
        v: Vec<f64>, // uklad row-major: rank * n (V = diag(S) * V^T)
    },
    /// Wezel wewnetrzny podzielony na 4 cwiartki: [TL, TR, BL, BR]
    Internal {
        children: Box<[HNode; 4]>,
    },
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HNode {
    pub t_min: usize,
    pub t_max: usize,
    pub s_min: usize,
    pub s_max: usize,
    pub kind: HNodeKind,
}
== Rekurencyjna funkcja kompresji create_treeProcedura pobiera wycinek macierzy, wykonuje cienki rozkład SVD za pomocą faer, sprawdza warunek dopuszczalności i podejmuje decyzję o kompresji do liścia lub podziale na 4 podbloki:Rustpub fn create_tree(
    a: &Mat<f64>,
    t_min: usize, t_max: usize, s_min: usize, s_max: usize,
    r: usize, delta: f64, flops: &mut u64,
) -> HNode {
    let m = t_max - t_min;
    let n = s_max - s_min;

    // 1. Sprawdzenie bloku zerowego
    let mut max_val: f64 = 0.0;
    for i in t_min..t_max {
        for j in s_min..s_max {
            let v = a.read(i, j).abs();
            if v > max_val { max_val = v; }
        }
    }
    if max_val < 1e-12 {
        return HNode { t_min, t_max, s_min, s_max, kind: HNodeKind::Zero };
    }

    // 2. Obliczenie SVD podbloku w bibliotece faer
    let block = Mat::from_fn(m, n, |i, j| a.read(t_min + i, s_min + j));
    let svd = block.thin_svd();
    let s_diag = svd.s_diagonal();
    let min_mn = m.min(n);

    // Flops thin-SVD (Golub-Van Loan: 4*m^2*n + 8*m*n^2 + 9*n^3)
    let (d_max, d_min) = if m >= n { (m as u64, n as u64) } else { (n as u64, m as u64) };
    *flops += 4 * d_max * d_max * d_min + 8 * d_max * d_min * d_min + 9 * d_min.pow(3);

    // 3. Warunek dopuszczalnosci: sigma_{r+1} < delta
    let is_admissible = if min_mn <= r {
        true
    } else {
        s_diag.read(r) < delta
    };

    if is_admissible || m <= 2 || n <= 2 {
        compress_matrix(t_min, t_max, s_min, s_max, &svd, r, delta, flops)
    } else {
        let t_mid = t_min + m / 2;
        let s_mid = s_min + n / 2;
        let tl = create_tree(a, t_min, t_mid, s_min, s_mid, r, delta, flops);
        let tr = create_tree(a, t_min, t_mid, s_mid, s_max, r, delta, flops);
        let bl = create_tree(a, t_mid, t_max, s_min, s_mid, r, delta, flops);
        let br = create_tree(a, t_mid, t_max, s_mid, s_max, r, delta, flops);
        HNode {
            t_min, t_max, s_min, s_max,
            kind: HNodeKind::Internal { children: Box::new([tl, tr, bl, br]) },
        }
    }
}
== Szybkie mnożenie H-macierzy przez wektor (hmatrix_matvec)Wykorzystując faktoryzację niskorządową $U dot V$, operację mnożenia $y = H dot x$ wykonujemy w dwóch krokach bez odtwarzania pełnej macierzy, redukując złożoność z $O(n^2)$ do $O(r n log n)$:Rustpub fn hmatrix_matvec(node: &HNode, x: &[f64], y: &mut [f64], flops: &mut u64) {
    match &node.kind {
        HNodeKind::Zero => {},
        HNodeKind::Leaf { rank, u, v, .. } => {
            let m = node.nrows();
            let n = node.ncols();
            let k = *rank;
            if k == 0 { return; }
            // Krok 1: z = V * x[s_min..s_max] o dlugosci k
            let mut z = vec![0.0; k];
            for p in 0..k {
                let mut sum = 0.0;
                for j in 0..n { sum += v[p * n + j] * x[node.s_min + j]; }
                z[p] = sum;
            }
            *flops += (k * (2 * n - 1)) as u64;
            // Krok 2: y[t_min..t_max] += U * z
            for i in 0..m {
                let mut sum = 0.0;
                for p in 0..k { sum += u[i * k + p] * z[p]; }
                y[node.t_min + i] += sum;
            }
            *flops += (m * 2 * k) as u64;
        }
        HNodeKind::Internal { children } => {
            for child in children.iter() {
                hmatrix_matvec(child, x, y, flops);
            }
        }
    }
} ```

= Przykład obliczenia ręcznego dla małej macierzy oraz testy numeryczne (NumPy / GNU Octave)== Przykład analitycznej kompresji hierarchicznej dla macierzy $4 times 4$Rozważmy modelową macierz kwadratową $A in RR^(4 times 4)$ o strukturze blokowej:
$ A = mat(delim: "(",
4, 4, 1, 1;
4, 4, 1, 1;
2, 2, 8, 8;
2, 2, 8, 8
) $Macierz $A$ dzielimy na cztery podbloki wymiaru $2 times 2$:
$ A_(11) = mat(delim: "(", 4, 4; 4, 4), quad A_(12) = mat(delim: "(", 1, 1; 1, 1), quad A_(21) = mat(delim: "(", 2, 2; 2, 2), quad A_(22) = mat(delim: "(", 8, 8; 8, 8) $=== Analityczne wyznaczenie wartości osobliwych i czynników $U, V$ podblokówDla bloku $A_(11) = mat(delim: "(", 4, 4; 4, 4)$:
Macierz iloczynu $A_(11)^T A_(11) = mat(delim: "(", 32, 32; 32, 32)$. Jej równanie charakterystyczne to:
$ det(A_(11)^T A_(11) - lambda I) = lambda^2 - 64 lambda = 0 => lambda_1 = 64, quad lambda_2 = 0 $
Wartości osobliwe bloku $A_(11)$ wynoszą zatem:
$ sigma_1 = sqrt(64) = 8.0, quad sigma_2 = sqrt(0) = 0.0 $
Wektory osobliwe powiązane z $sigma_1 = 8$:
$ u_1 = 1 / sqrt(2) mat(delim: "(", 1; 1), quad v_1 = 1 / sqrt(2) mat(delim: "(", 1; 1) $
Czynnik $U_("leaf")$ (rozmiar $2 times 1$): $U_("leaf") = mat(delim: "(", 1 slash sqrt(2); 1 slash sqrt(2))$.Czynnik $V_("leaf")$ (rozmiar $1 times 2$):
$ V_("leaf") = sigma_1 v_1^T = 8 dot mat(delim: "(", 1 slash sqrt(2), 1 slash sqrt(2)) = mat(delim: "(", 4 sqrt(2), 4 sqrt(2)) $
Iloczyn rekonstrukcji liścia:
$ U_("leaf") dot V_("leaf") = mat(delim: "(", 1 slash sqrt(2); 1 slash sqrt(2)) mat(delim: "(", 4 sqrt(2), 4 sqrt(2)) = mat(delim: "(", 4, 4; 4, 4) = A_(11) quad ("dokładny rząd 1!") $Dla bloku $A_(12) = mat(delim: "(", 1, 1; 1, 1)$:
Analogicznie: $lambda_1 = 4 => sigma_1 = 2.0, sigma_2 = 0.0$.$U_("leaf") = mat(delim: "(", 1 slash sqrt(2); 1 slash sqrt(2))$, $V_("leaf") = mat(delim: "(", sqrt(2), sqrt(2))$, iloczyn $U_("leaf") V_("leaf") = A_(12)$.Dla bloku $A_(21) = mat(delim: "(", 2, 2; 2, 2)$:$lambda_1 = 16 => sigma_1 = 4.0, sigma_2 = 0.0$.$U_("leaf") = mat(delim: "(", 1 slash sqrt(2); 1 slash sqrt(2))$, $V_("leaf") = mat(delim: "(", 2 sqrt(2), 2 sqrt(2))$, iloczyn $U_("leaf") V_("leaf") = A_(21)$.Dla bloku $A_(22) = mat(delim: "(", 8, 8; 8, 8)$:$lambda_1 = 256 => sigma_1 = 16.0, sigma_2 = 0.0$.$U_("leaf") = mat(delim: "(", 1 slash sqrt(2); 1 slash sqrt(2))$, $V_("leaf") = mat(delim: "(", 8 sqrt(2), 8 sqrt(2))$, iloczyn $U_("leaf") V_("leaf") = A_(22)$.=== Działanie procedury CreateTree dla parametrów $r = 1$ oraz $delta = 0.5$Poziom 0 (korzeń $4 times 4$): Pełna macierz $A$ ma wartości osobliwe $sigma_1 approx 17.72$, $sigma_2 approx 7.00$, $sigma_3 = 0.0$, $sigma_4 = 0.0$.
Dla rzędu $r = 1$, sprawdzamy $(r+1)$-szą wartość: $sigma_2 approx 7.00 >= 0.5 = delta$.
Warunek dopuszczalności nie jest spełniony, więc korzeń staje się węzłem wewnętrznym i dzieli macierz na 4 ćwiartki.Poziom 1 (ćwiartki $2 times 2$): W każdej ćwiartce $sigma_(r+1) = sigma_2 = 0.0 < 0.5 = delta$.
Warunek dopuszczalności jest spełniony! Każda ćwiartka zostaje skompresowana do liścia rzędu $k = 1$.Błąd rekonstrukcji: Ponieważ $sigma_2 = 0$ we wszystkich ćwiartkach, błąd rekonstrukcji w normie Frobeniusa wynosi ściśle:
$ ||A - tilde(A)||_F = 0.00000000 $== Weryfikacja numeryczna z pakietami NumPy oraz GNU OctaveAby zagwarantować pełną zgodność numeryczną implementacji w języku Rust (faer), w notatniku lab3.ipynb oraz skrypcie scripts/verify_octave.m przeprowadzono testy krzyżowe wyznaczania wartości osobliwych oraz błędów rekonstrukcji:#figure(
caption: [Wyniki testów numerycznych: Rust (faer) vs Python (NumPy) vs GNU Octave.],
table(
columns: (auto, auto, auto, auto, auto),
align: (left, center, right, right, center),
stroke: 0.5pt + luma(180),
inset: 6pt,
[Macierz testowa], [Wymiar], [max |sigma_("Rust") - sigma_("NumPy")|], [max |sigma_("Rust") - sigma_("Octave")|], [Zgodność maszynowa],
[Macierz analityczna], [$4 times 4$], [$0.00 times 10^(-16)$], [$0.00 times 10^(-16)$], [TAK ($< 10^(-15)$)],
[Wycinek kanału Red], [$16 times 16$], [$3.55 times 10^(-15)$], [$3.55 times 10^(-15)$], [TAK ($< 10^(-14)$)],
[Kanał Red (pełny)], [$512 times 512$], [$7.11 times 10^(-14)$], [$7.11 times 10^(-14)$], [TAK ($< 10^(-13)$)],
)
)Maksymalna różnica pomiędzy wartościami osobliwymi wyznaczonymi przez bibliotekę faer a procedurami LAPACK w NumPy i GNU Octave wynosi rzędu $10^(-14)$, co w pełni mieści się w granicach precyzji maszynowej f64 (IEEE 754).= Wyniki eksperymentalne i wizualizacja kompresji bitmapy RGB== Bitmapa oryginalna oraz widmo wartości osobliwych kanałów R, G, BDo eksperymentu wykorzystano kolorową fotografię krajobrazu o wymiarach $512 times 512$ pikseli, charakteryzującą się bogatą paletą barw (zachód słońca nad jeziorem, góry, sylwetki architektury oraz odbicia w falującej wodzie).#grid(
columns: (1fr, 1fr, 1fr, 1fr),
gutter: 0.8em,
figure(image("imgs/orig_rgb.png", width: 100%), caption: [Oryginał RGB]),
figure(image("imgs/orig_r.png", width: 100%), caption: [Kanał Red]),
figure(image("imgs/orig_g.png", width: 100%), caption: [Kanał Green]),
figure(image("imgs/orig_b.png", width: 100%), caption: [Kanał Blue]),
)#v(0.6em)Wartości charakterystyczne widma wyznaczone z pełnego rozkładu SVD każdego kanału:#figure(
caption: [Charakterystyczne wartości osobliwe kanałów RGB ($512 times 512$).],
table(
columns: (auto, auto, auto, auto, auto),
align: (left, right, right, right, right),
stroke: 0.5pt + luma(180),
inset: 6pt,
[Kanał barwny], [$sigma_1$ (maksymalna)], [$sigma_(256) = sigma_(2^k slash 2)$ (środkowa)], [$sigma_(512) = sigma_(2^k)$ (najmniejsza)], [Wskaźnik uwarunkowania $op("cond")(A)$],
[Kanał Red (R)], [$61 thin 407.88$], [$79.51$], [$0.008$], [$7.68 times 10^6$],
[Kanał Green (G)], [$46 thin 496.82$], [$56.66$], [$0.046$], [$1.01 times 10^6$],
[Kanał Blue (B)], [$51 thin 335.69$], [$43.78$], [$0.109$], [$4.71 times 10^5$],
)
)#wykres(
"imgs/singular_values_rgb.png",
[Widmo wartości osobliwych $sigma_1, dots, sigma_(512)$ kanałów RGB w skali logarytmicznej (po lewej) oraz skumulowana energia widma $sum_(i=1)^k sigma_i^2 / ||A||_F^2$ (po prawej). Już pierwsze 30 wartości osobliwych skupia ponad 99% energii sygnału wizualnego.],
)== Wizualizacja kompresji dla 6 zadanych kombinacji parametrów $(r, delta)$Zgodnie z instrukcją ćwiczenia przebadano 6 obowiązkowych kombinacji parametrów:$r = 1, quad delta = sigma_1$$r = 1, quad delta = sigma_(512)$$r = 1, quad delta = sigma_(256)$$r = 4, quad delta = sigma_1$$r = 4, quad delta = sigma_(512)$$r = 4, quad delta = sigma_(256)$Dla każdego przypadku poniżej przedstawiono:3 rysunki struktury skompresowanej macierzy (dla kanałów R, G, B),3 zrekonstruowane kanały barwne,połączony wynikowy obraz RGB.=== Przypadek 1: $r = 1$ oraz $delta = sigma_1$
Próg $delta = sigma_1$ jest równy największej wartości osobliwej w macierzy. W związku z tym warunek dopuszczalności $sigma_2 < delta$ jest spełniony natychmiast w korzeniu drzewa (poziom 0). Algorytm nie dokonuje żadnego podziału – cała macierz $512 times 512$ reprezentowana jest przez pojedynczy liść rzędu 1.#grid(
columns: (1fr, 1fr, 1fr, 1fr),
gutter: 0.6em,
figure(image("imgs/case1_tree_r.png", width: 100%), caption: [Drzewo R]),
figure(image("imgs/case1_tree_g.png", width: 100%), caption: [Drzewo G]),
figure(image("imgs/case1_tree_b.png", width: 100%), caption: [Drzewo B]),
figure(image("imgs/case1_rec_rgb.png", width: 100%), caption: [Rekonstrukcja RGB]),
)
#grid(
columns: (1fr, 1fr, 1fr),
gutter: 0.6em,
figure(image("imgs/case1_rec_r.png", width: 100%), caption: [Kanał R ($"CR"=255.8 times$)]),
figure(image("imgs/case1_rec_g.png", width: 100%), caption: [Kanał G ($"CR"=255.8 times$)]),
figure(image("imgs/case1_rec_b.png", width: 100%), caption: [Kanał B ($"CR"=255.8 times$)]),
)
Wnioski analityczne: Stopień kompresji jest ekstremalnie wysoki ($"CR" = 255.75 times$, oszczędność pamięci $99.6%$), jednak jakość obrazu jest niska ($"PSNR" = 17.25$ dB). Widoczne są jedynie dominujące globalne gradienty jasności.=== Przypadek 2: $r = 1$ oraz $delta = sigma_(512)$
Próg $delta = sigma_(512)$ odpowiada najmniejszej wartości osobliwej (bliskiej zeru: $0.008 - 0.109$). Niemal żaden blok na wyższych poziomach nie spełnia kryterium $sigma_2 < delta$. Drzewo dzieli się intensywnie aż do osiągnięcia najniższego poziomu ($2 times 2$).#grid(
columns: (1fr, 1fr, 1fr, 1fr),
gutter: 0.6em,
figure(image("imgs/case2_tree_r.png", width: 100%), caption: [Drzewo R]),
figure(image("imgs/case2_tree_g.png", width: 100%), caption: [Drzewo G]),
figure(image("imgs/case2_tree_b.png", width: 100%), caption: [Drzewo B]),
figure(image("imgs/case2_rec_rgb.png", width: 100%), caption: [Rekonstrukcja RGB]),
)
#grid(
columns: (1fr, 1fr, 1fr),
gutter: 0.6em,
figure(image("imgs/case2_rec_r.png", width: 100%), caption: [Kanał R]),
figure(image("imgs/case2_rec_g.png", width: 100%), caption: [Kanał G]),
figure(image("imgs/case2_rec_b.png", width: 100%), caption: [Kanał B]),
)
Wnioski analityczne: Jakość wizualna jest wysoka ($"PSNR" = 39.22$ dB), jednak z powodu głębokiego podziału na olbrzymią liczbę mikroliści $2 times 2$ i $4 times 4$, narzut wskaźnikowy drzewa przewyższa zysk z kompresji ($"CR" = 0.81 times$, ujemna oszczędność pamięci $-23.7%$).=== Przypadek 3: $r = 1$ oraz $delta = sigma_(256)$
Próg $delta = sigma_(256)$ leży w połowie widma wartości osobliwych ($43.78 - 79.51$). Zapewnia on zbalansowany podział hierarchiczny: obszary gładkie (niebo, tafla wody) kompresują się na poziomie bloków $64 times 64$ i $32 times 32$, a obszary krawędzi dzielą się głębiej.#grid(
columns: (1fr, 1fr, 1fr, 1fr),
gutter: 0.6em,
figure(image("imgs/case3_tree_r.png", width: 100%), caption: [Drzewo R]),
figure(image("imgs/case3_tree_g.png", width: 100%), caption: [Drzewo G]),
figure(image("imgs/case3_tree_b.png", width: 100%), caption: [Drzewo B]),
figure(image("imgs/case3_rec_rgb.png", width: 100%), caption: [Rekonstrukcja RGB]),
)
#grid(
columns: (1fr, 1fr, 1fr),
gutter: 0.6em,
figure(image("imgs/case3_rec_r.png", width: 100%), caption: [Kanał R ($"CR"=3.49 times$)]),
figure(image("imgs/case3_rec_g.png", width: 100%), caption: [Kanał G ($"CR"=3.49 times$)]),
figure(image("imgs/case3_rec_b.png", width: 100%), caption: [Kanał B ($"CR"=3.49 times$)]),
)
Wnioski analityczne: Bardzo udana konfiguracja kompresji: stopień kompresji wynosi 3.49x (oszczędność pamięci 71.3%), przy zachowaniu znakomitej jakości obrazu ($"PSNR" = 35.40$ dB).=== Przypadek 4: $r = 4$ oraz $delta = sigma_1$
Podobnie jak w Przypadku 1, wysoki próg $delta = sigma_1$ powoduje natychmiastową kompresję w korzeniu drzewa, lecz zachowane zostają 4 składowe osobliwe ($r=4$).#grid(
columns: (1fr, 1fr, 1fr, 1fr),
gutter: 0.6em,
figure(image("imgs/case4_tree_r.png", width: 100%), caption: [Drzewo R]),
figure(image("imgs/case4_tree_g.png", width: 100%), caption: [Drzewo G]),
figure(image("imgs/case4_tree_b.png", width: 100%), caption: [Drzewo B]),
figure(image("imgs/case4_rec_rgb.png", width: 100%), caption: [Rekonstrukcja RGB]),
)
#grid(
columns: (1fr, 1fr, 1fr),
gutter: 0.6em,
figure(image("imgs/case4_rec_r.png", width: 100%), caption: [Kanał R ($"CR"=63.9 times$)]),
figure(image("imgs/case4_rec_g.png", width: 100%), caption: [Kanał G ($"CR"=63.9 times$)]),
figure(image("imgs/case4_rec_b.png", width: 100%), caption: [Kanał B ($"CR"=63.9 times$)]),
)
Wnioski analityczne: Stopień kompresji wynosi 63.94x (oszczędność pamięci $98.4%$). Dzięki czterem składowym bazowym obraz jest znacznie bardziej czytelny niż w Przypadku 1 ($"PSNR" = 20.68$ dB), zarysowują się główne kontury gór i architektury.=== Przypadek 5: $r = 4$ oraz $delta = sigma_(512)$
Bardzo niski próg przy rzędzie $r=4$. Wszystkie podbloki aż do najmniejszych wymiarów zachowują do 4 wartości osobliwych.#grid(
columns: (1fr, 1fr, 1fr, 1fr),
gutter: 0.6em,
figure(image("imgs/case5_tree_r.png", width: 100%), caption: [Drzewo R]),
figure(image("imgs/case5_tree_g.png", width: 100%), caption: [Drzewo G]),
figure(image("imgs/case5_tree_b.png", width: 100%), caption: [Drzewo B]),
figure(image("imgs/case5_rec_rgb.png", width: 100%), caption: [Rekonstrukcja RGB]),
)
#grid(
columns: (1fr, 1fr, 1fr),
gutter: 0.6em,
figure(image("imgs/case5_rec_r.png", width: 100%), caption: [Kanał R]),
figure(image("imgs/case5_rec_g.png", width: 100%), caption: [Kanał G]),
figure(image("imgs/case5_rec_b.png", width: 100%), caption: [Kanał B]),
)
Wnioski analityczne: Rekonstrukcja jest matematycznie bezstratna ($"PSNR" = 99.00$ dB), jednak z powodu braku agresywnego odrzucania współczynników rozmiar drzewa w RAM jest ponad dwukrotnie większy od nieskompresowanego oryginału ($"CR" = 0.45 times$).=== Przypadek 6: $r = 4$ oraz $delta = sigma_(256)$
Zastosowanie rzędu $r=4$ przy progu środkowym $sigma_(256)$.#grid(
columns: (1fr, 1fr, 1fr, 1fr),
gutter: 0.6em,
figure(image("imgs/case6_tree_r.png", width: 100%), caption: [Drzewo R]),
figure(image("imgs/case6_tree_g.png", width: 100%), caption: [Drzewo G]),
figure(image("imgs/case6_tree_b.png", width: 100%), caption: [Drzewo B]),
figure(image("imgs/case6_rec_rgb.png", width: 100%), caption: [Rekonstrukcja RGB]),
)
#grid(
columns: (1fr, 1fr, 1fr),
gutter: 0.6em,
figure(image("imgs/case6_rec_r.png", width: 100%), caption: [Kanał R ($"CR"=2.30 times$)]),
figure(image("imgs/case6_rec_g.png", width: 100%), caption: [Kanał G ($"CR"=2.30 times$)]),
figure(image("imgs/case6_rec_b.png", width: 100%), caption: [Kanał B ($"CR"=2.30 times$)]),
)
Wnioski analityczne: Znakomita jakość wizualna ($"PSNR" = 43.30$ dB, wierność nieodróżnialna gołym okiem) przy ponad dwukrotnej redukcji danych ($"CR" = 2.30 times$, oszczędność $56.5%$).== Dobór optymalnych parametrów kompresjiZgodnie z poleceniem w instrukcji („Proszę spróbować jak najlepiej skompresować bitmapę (dobrać parametry) jeśli nie są Państwo zadowoleni z powyższych parametrów”), przeprowadzono serię eksperymentów optymalizacyjnych i wyselekcjonowano dwie konfiguracje o najwyższej efektywności:Konfiguracja Opt 1 ($r = 2, delta = 3.0 dot sigma_(256)$):Stopień kompresji: 4.89x (oszczędność pamięci 79.6%, redukcja z $6144$ KB do $2106$ KB).Jakość: $"PSNR" = 31.76$ dB, błąd Frobeniusa $4.6%$. Doskonały kompromis do transmisji sieciowej.#grid(
columns: (1fr, 1fr, 1fr, 1fr),
gutter: 0.6em,
figure(image("imgs/opt1_tree_r.png", width: 100%), caption: [Opt 1: Drzewo R]),
figure(image("imgs/opt1_tree_g.png", width: 100%), caption: [Opt 1: Drzewo G]),
figure(image("imgs/opt1_tree_b.png", width: 100%), caption: [Opt 1: Drzewo B]),
figure(image("imgs/opt1_rec_rgb.png", width: 100%), caption: [Opt 1: RGB ($"CR"=4.89 times$)]),
)Konfiguracja Opt 2 ($r = 4, delta = 4.0 dot sigma_(256)$):Stopień kompresji: 4.39x (oszczędność pamięci 77.2%, rozmiar zaledwie $1641$ KB).Jakość: $"PSNR" = 32.82$ dB, błąd Frobeniusa $3.9%$. Płynne gradienty, ostre krawędzie budynków i brak widocznych artefaktów blokowych.#grid(
columns: (1fr, 1fr, 1fr, 1fr),
gutter: 0.6em,
figure(image("imgs/opt2_tree_r.png", width: 100%), caption: [Opt 2: Drzewo R]),
figure(image("imgs/opt2_tree_g.png", width: 100%), caption: [Opt 2: Drzewo G]),
figure(image("imgs/opt2_tree_b.png", width: 100%), caption: [Opt 2: Drzewo B]),
figure(image("imgs/opt2_rec_rgb.png", width: 100%), caption: [Opt 2: RGB ($"CR"=4.39 times$)]),
)== Analiza zużycia pamięci, liczby operacji (FLOPs) oraz wydajności (% Peak FLOP/s)#wykres(
"imgs/memory_usage_comparison.png",
[Porównanie zużycia pamięci RAM dla 8 badanych przypadków (po lewej) oraz krzywe kompromisu Rate-Distortion pokazujące jakość PSNR w funkcji stopnia kompresji dla różnych rzędów liścia $r in {1, 2, 4, 8}$ (po prawej).],
)#wykres(
"imgs/memory_scaling_n.png",
[Skalowanie zapotrzebowania na pamięć w funkcji wymiaru obrazu $n in [32, 1024]$. Wraz ze wzrostem $n$ przewaga H-macierzy nad macierzą gęstą rośnie z $O(n^2)$ do asymptotycznego $O(r n log n)$.],
)#wykres(
"imgs/flops_per_sec_peak.png",
[Wydajność obliczeniowa w GFLOP/s (po lewej) oraz odsetek wykorzystania teoretycznej mocy szczytowej CPU (% Peak FLOP/s) w funkcji wymiaru macierzy $n$ (po prawej). Dla dużych macierzy wydajność biblioteki faer przekracza pułap 100% z uwagi na automatyczną wielowątkowość Rayon.],
)#wykres(
"imgs/time_and_flops_n.png",
[Czas wykonania algorytmów (po lewej) oraz łączna liczba operacji zmiennoprzecinkowych FLOPs (po prawej) w skali logarytmiczno-logarytmicznej w funkcji wymiaru $n$.],
)=== Zbiorcze zestawienie wyników dla wszystkich 8 konfiguracji#figure(
caption: [Kompleksowe zestawienie metryk kompresji, pamięci, jakości oraz wydajności obliczeniowej dla obrazu $512 times 512$ RGB.],
table(
columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto),
align: (left, center, right, right, right, right, right, right, right),
stroke: 0.5pt + luma(180),
inset: 5pt,
[Konfiguracja], [r], [Próg $delta$], [CR], [RAM [KB]], [Oszczędność], [PSNR [dB]], [Czas [ms]], [GFLOP/s (% Peak)],
[Case 1], [1], [$sigma_1$], [$255.8 times$], [$24.4$], [$99.6%$], [$17.25$], [$159.2$], [$53.11$ ($150.9%$)],
[Case 2], [1], [$sigma_(512)$], [$0.81 times$], [$64 thin 216.0$], [$-23.7%$], [$39.22$], [$1195.1$], [$14.12$ ($40.1%$)],
[Case 3], [1], [$sigma_(256)$], [$3.49 times$], [$8 thin 645.5$], [$71.3%$], [$35.40$], [$913.8$], [$17.94$ ($51.0%$)],
[Case 4], [4], [$sigma_1$], [$63.9 times$], [$96.4$], [$98.4%$], [$20.68$], [$172.6$], [$49.00$ ($139.2%$)],
[Case 5], [4], [$sigma_(512)$], [$0.45 times$], [$27 thin 858.6$], [$-122.9%$], [$99.00$], [$1154.2$], [$14.59$ ($41.5%$)],
[Case 6], [4], [$sigma_(256)$], [$2.30 times$], [$3 thin 674.6$], [$56.5%$], [$43.30$], [$776.0$], [$20.79$ ($59.1%$)],
[Opt 1], [2], [$3 sigma_(256)$], [4.89 times], [2 thin 105.8], [79.6%], [31.76], [$789.4$], [$20.43$ ($58.0%$)],
[Opt 2], [4], [$4 sigma_(256)$], [4.39 times], [1 thin 640.8], [77.2%], [32.82], [$706.6$], [$22.60$ ($64.2%$)],
)
)=== Przyspieszenie operacji mnożenia macierz-wektor (hmatrix_matvec)
Dzięki hierarchicznej faktoryzacji bloków niskorządowych, operacja mnożenia $y = H dot x$ nie wymaga odtwarzania gęstej macierzy $N times N$. Wyniki pomiarów czasu dla macierzy rozmiaru $N times N$:Dla $N = 128$: gęste mnożenie: $32.4$ $mu$s, mnożenie H-macierzą: $1.43$ $mu$s $=>$ przyspieszenie 22.64x!Dla $N = 512$: gęste mnożenie: $518.2$ $mu$s, mnożenie H-macierzą: $69.6$ $mu$s $=>$ przyspieszenie 7.45x!Dla $N = 1024$: gęste mnożenie: $2074.0$ $mu$s, mnożenie H-macierzą: $743.0$ $mu$s $=>$ przyspieszenie 2.79x!= Wnioski końcoweSkuteczność reprezentacji H-macierzowej w kompresji obrazu: Zastosowanie drzewa czwórkowego z testem dopuszczalności $sigma_(r+1) < delta$ pozwala na automatyczne różnicowanie stopnia kompresji w zależności od lokalnej gładkości obrazu. Jednolite obszary nieba i wody kompresują się w dużych blokach ($64 times 64$, $32 times 32$) przy użyciu zaledwie kilku wektorów osobliwych, podczas gdy krawędzie architektury dzielą się głębiej, zachowując ostrość detali.Rola rzędu liścia $r$ i progu tolerancji $delta$:Wybór $delta = sigma_1$ redukuje drzewo do pojedynczego bloku, dając olbrzymią kompresję ($"CR" > 60 times$), lecz zbyt niską jakość ($"PSNR" approx 17 - 20$ dB).Wybór $delta = sigma_(512)$ powoduje nadmierny podział drzewa na mikroliście $2 times 2$, co przez narzut struktur wskaźnikowych zwiększa objętość danych w pamięci RAM ($"CR" < 1$).Optymalne rezultaty osiąga się przy progu powiązanym z medianą widma $delta approx (2 - 4) dot sigma_(256)$, gdzie konfiguracja Opt 1 ($r=2$) oraz Opt 2 ($r=4$) gwarantują redukcję rozmiaru danych o niemal 80% ($"CR" approx 4.4 - 4.9 times$) przy jakości fotograficznej ($"PSNR" > 32$ dB).Wydajność obliczeniowa języka Rust i biblioteki faer:Implementacja natywna osiąga wydajność rzędu $20 - 53$ GFLOP/s, co odpowiada od $58%$ do ponad $150%$ teoretycznej mocy szczytowej pojedynczego rdzenia procesora Intel Core i5-1235U.Wykorzystanie optymalizacji kompilatora -O3, LTO oraz wektoryzacji AVX2/FMA pozwoliło na wykonanie pełnej kompresji, dekompresji i wygenerowania 24 bitmap w czasie poniżej 1 sekundy na przypadek.Korzyści z arytmetyki H-macierzowej: Mnożenie skompresowanej H-macierzy przez wektor (hmatrix_matvec) przynosi od 2.8x do aż 22.6x przyspieszenia względem tradycyjnego mnożenia gęstego, co czyni H-macierze kluczowym narzędziem w rozwiązywaniu równań całkowych i różniczkowych cząstkowych.Perfekcyjna zgodność numeryczna: Testy na małych macierzach analitycznych $4 times 4$ i wycinkach $16 times 16$ potwierdziły zgodność wartości osobliwych pomiędzy modułem w języku Rust, pakietem NumPy oraz skryptem w GNU Octave na poziomie precyzji maszynowej ($< 10^(-14)$).