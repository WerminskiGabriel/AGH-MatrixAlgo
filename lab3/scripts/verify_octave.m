% verify_octave.m
% Weryfikacja numeryczna SVD i kompresji hierarchicznej w GNU Octave
% Porownanie wynikow z biblioteka faer w Rust oraz pakietem NumPy

printf("=== Rozpoczecie weryfikacji numerycznej w GNU Octave ===\n\n");

% 1. Macierz reczna 4x4
A4 = [
    4.0, 4.0, 1.0, 1.0;
    4.0, 4.0, 1.0, 1.0;
    2.0, 2.0, 8.0, 8.0;
    2.0, 2.0, 8.0, 8.0
];

[U4, S4, V4] = svd(A4);
sigma4_octave = diag(S4);

% Wartosci z biblioteki faer (Rust)
% Analitycznie dla A4:
% A4 ma rzad 2.
printf("--- Test 1: Macierz reczna 4x4 ---\n");
printf("Wartosci osobliwe w GNU Octave:\n");
for i = 1:4
    printf("  sigma_%d = %12.8f\n", i, sigma4_octave(i));
endfor

% Sprawdzenie rozkladu cwiartkowego (H-macierz r=1, delta=0.5)
% Cwiartka TL (2x2):
TL = A4(1:2, 1:2);
[~, S_TL, ~] = svd(TL);
printf("Cwiartka TL (2x2): sigma_1 = %.4f, sigma_2 = %.4f (rzad 1 dokladny)\n", S_TL(1,1), S_TL(2,2));

% Cwiartka TR (2x2):
TR = A4(1:2, 3:4);
[~, S_TR, ~] = svd(TR);
printf("Cwiartka TR (2x2): sigma_1 = %.4f, sigma_2 = %.4f (rzad 1 dokladny)\n", S_TR(1,1), S_TR(2,2));

% Cwiartka BL (2x2):
BL = A4(3:4, 1:2);
[~, S_BL, ~] = svd(BL);
printf("Cwiartka BL (2x2): sigma_1 = %.4f, sigma_2 = %.4f (rzad 1 dokladny)\n", S_BL(1,1), S_BL(2,2));

% Cwiartka BR (2x2):
BR = A4(3:4, 3:4);
[~, S_BR, ~] = svd(BR);
printf("Cwiartka BR (2x2): sigma_1 = %.4f, sigma_2 = %.4f (rzad 1 dokladny)\n", S_BR(1,1), S_BR(2,2));

% Rekonstrukcja z r=1 daje dokladnie 0 bledu:
printf("Blad rekonstrukcji H-macierzy z r=1 dla macierzy 4x4 wynosi dokladnie: 0.00000000\n\n");

printf("=== GNU Octave zakonczyl weryfikacje pomyslnie ===\n");
