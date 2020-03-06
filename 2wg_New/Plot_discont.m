
F = linspace(6e9, 12e9, 35);

%% standard mesh


standard = read(rfdata.data,'2wg_standard_mesh.s20p');
S_std = extract(standard,'S_PARAMETERS');


%% Fine mesh 

fine = read(rfdata.data,'2wg_fine_mesh.s20p');
S_fine = extract(fine,'S_PARAMETERS');


%% Super Fine mesh

super_fine = read(rfdata.data,'2wg_super_fine_mesh.s20p');
S_super_fine = extract(super_fine,'S_PARAMETERS');

%% Plots

figure;

plot(F * 1e-9, db(abs(squeeze(S_std(11, 15, :))))/2, 'LineWidth', 1); grid on;
hold on;

plot(F * 1e-9, db(abs(squeeze(S_fine(11, 15, :))))/2, 'LineWidth', 1); grid on;
hold on;

plot(F * 1e-9, db(abs(squeeze(S_super_fine(11, 15, :))))/2, 'LineWidth', 1); grid on;


% xlim([9.46 9.48]);
% ylim([-21 -20.5]);
xlim([6 10.5])