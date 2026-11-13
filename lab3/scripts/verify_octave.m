% verify_octave.m
% Weryfikacja wartosci osobliwych: Rust (faer) vs GNU Octave.
% Uruchamiac z katalogu lab3:  octave --no-gui --quiet scripts/verify_octave.m
%
% Dane z Rusta (results/small_test.json, results/results.json) sa porownywane
% z wynikiem svd() z Octave (LAPACK) dla tych samych macierzy.

printf("=== Weryfikacja numeryczna SVD: Rust (faer) vs GNU Octave ===\n\n");

small = jsondecode(fileread("results/small_test.json"));
res   = jsondecode(fileread("results/results.json"));
img   = double(imread("assets/input_512.png"));

names  = {};
dims   = {};
absd   = [];
reld   = [];

% --- macierze testowe (4x4 oraz wycinek 16x16) ---
if iscell(small)
  tests = small;
else
  tests = num2cell(small);
end
for k = 1:numel(tests)
  t = tests{k};
  A = double(t.matrix);
  s_rust = double(t.singular_values(:));
  s_oct  = svd(A);
  names{end+1} = t.name;
  dims{end+1}  = sprintf("%dx%d", rows(A), columns(A));
  absd(end+1)  = max(abs(s_rust - s_oct));
  reld(end+1)  = absd(end) / s_oct(1);
end

% --- pelne kanaly R, G, B (512x512) ---
chan = {"r", "g", "b"};
for c = 1:3
  A = img(:, :, c);
  s_rust = double(res.(["singular_values_" chan{c}]).values(:));
  s_oct  = svd(A);
  names{end+1} = ["channel_" chan{c}];
  dims{end+1}  = sprintf("%dx%d", rows(A), columns(A));
  absd(end+1)  = max(abs(s_rust - s_oct));
  reld(end+1)  = absd(end) / s_oct(1);
end

printf("%-20s %-9s %-28s %-26s\n", "Macierz", "Wymiar", "max|sigma_Rust-sigma_Octave|", "blad wzgledny / sigma_1");
for k = 1:numel(names)
  printf("%-20s %-9s %-28.2e %-26.2e\n", names{k}, dims{k}, absd(k), reld(k));
end

% zapis do CSV
fid = fopen("results/octave_compare.csv", "w");
fprintf(fid, "matrix,dim,max_abs_diff,rel_diff_sigma1\n");
for k = 1:numel(names)
  fprintf(fid, "%s,%s,%.6e,%.6e\n", names{k}, dims{k}, absd(k), reld(k));
end
fclose(fid);

printf("\nZapisano: results/octave_compare.csv\n");
printf("=== Koniec ===\n");
