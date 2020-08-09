c0 = 3e8;

R = [2e-2 4e-2];
Len = 5e-2;

F = linspace(5e9, 21e9, 35);

[STT, STR, SRT, SRR, num] = GSM_N_Vfinal(R, Len, F, 20);

figure;
plot(F .* 1e-9, db(abs(STT(:, 7, 7))), 'LineWidth', 2); grid on;
figure;
plot(F .* 1e-9, db(abs(STR(:, 7, 7))), 'LineWidth', 2); grid on;

figure;
plot(F .* 1e-9, angle(STT(:, 7, 7)), 'LineWidth', 2); grid on;
figure;
plot(F .* 1e-9, angle(STR(:, 7, 7)), 'LineWidth', 2); grid on;
