clear;
F = 5e9;
c0 = 3e8;


er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability



R = 2e-2;

n = 25:5:50;

for i = 1:length(n)

R_FF = n(i) .* R;

r = [R R_FF];

er = [1 1-10j];
mur = [1 1];

F = 5e9;

lamb = c0/F;

for k = 1:length(r)
    f =  fc(r(k), er(k), mur(k));
    N_i  =  find(f < F);
    N_(k) = length(N_i);
end


for m = 1:length(r)-1
    X_(i).X_til = Inner_p2(1:1:N_(m), 1:1:N_(m+1), r(m+1), r(m), er(m+1), mur(m+1), er(m), mur(m));
    [STT(i).stt, STR(i).str, SRT(i).srt, SRR(i).srr] = GSM_V2(1:1:N_(m), 1:1:N_(m+1), F, r(m+1), r(m), er(m+1), mur(m+1), er(m), mur(m), (X_(i).X_til));
    Gamma_TE11(i) = STT(i).stt(1, 1);
end
end

figure;
plot(n, abs(Gamma_TE11), 'LineWidth', 2);
grid on;

figure;
plot(n, angle(Gamma_TE11), 'LineWidth', 2);
grid on;