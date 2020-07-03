clear L
out = load('ga_allvar_ms3serv2.mat');

out = out.fmin2;

c0 = 3e8;
F = 5e9;
lamb = c0./F;
N = 16;

l1 = lamb/4;

R = out.r(1:N);

for i = N+1:length(out.r)
    L(i-N) = sum(out.r(N+1:i)) + l1;
end
% L(1) = l1;
L = [l1 L];

figure; plot(L, R, 'Color', [0.6350, 0.0780, 0.1840] , 'LineWidth', 2); grid on; hold on; plot(L, -R, 'Color', [0.6350, 0.0780, 0.1840] , 'LineWidth', 2); 