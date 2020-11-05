clear;
c0 = 3e8;
F = 5e9;
lamb = c0./F;
N = 5;
data = load('ps_V2L_ms3serv2_fl_V2.mat');
fmin = data.fmin2;

R = [fmin.r(1) fmin.r(1:N).'];

% R = [0.02283 0.02283 0.06313 0.07772 0.09971 0.1467]; % From GA algo


% R = [0.02 0.02 0.02667 0.03333 0.04 0.04667 0.05333 0.06 0.06667 0.07333 0.09562 0.08667 0.09333 0.1009 0.1067 0.1156...
%     0.1491];
% N = length(R) - 1;
% Len = 0.3912;
% Len = fmin.r(N + 1);

Len = fmin.r(N+1);

l1 = lamb./4;
L = [l1 ones(1, N - 1)./N .* Len];

for i = 1:N
    L_axis(i) = sum(L(1:i));
end
L_axis = [0 L_axis];



c0 = 3e8;
F = 5e9;
lamb = c0./F;
N = 5;
data = load('ga_V2L_ms3serv2_fl_V2.mat');
fmin = data.fmin2;

R2 = [fmin.r(1) fmin.r(1:N)];

% R = [0.02283 0.02283 0.06313 0.07772 0.09971 0.1467]; % From GA algo


% R = [0.02 0.02 0.02667 0.03333 0.04 0.04667 0.05333 0.06 0.06667 0.07333 0.09562 0.08667 0.09333 0.1009 0.1067 0.1156...
%     0.1491];
% N = length(R) - 1;
% Len = 0.3912;
% Len = fmin.r(N + 1);

Len2 = fmin.r(N+1);

l1 = lamb./4;
L = [l1 ones(1, N - 1)./N .* Len2];

for i = 1:N
    L_axis2(i) = sum(L(1:i));
end
L_axis2 = [0 L_axis2];

figure;
plot(L_axis, R, 'LineWidth', 2, 'Color', [0.6350, 0.0780, 0.1840] );

hold on;

plot(L_axis2, R2, 'LineWidth', 2, 'Color', [0.1840, 0.0780, 0.6350]);
hold on;

plot(L_axis, -R, 'LineWidth', 2, 'Color', [0.6350, 0.0780, 0.1840] );

hold on;
plot(L_axis2, -R2, 'LineWidth', 2, 'Color', [0.1840, 0.0780, 0.6350] );


legend({['Fmincon R_{ap} = ', num2str(R(end)./lamb), ' \lambda ', 'L = ', num2str(Len./lamb), ' \lambda'], ...
    ['GA R_{ap} = ', num2str(R2(end)./lamb), ' \lambda ', 'L = ', num2str(Len2./lamb), ' \lambda']},...
    'FontSize', 12, 'FontWeight', 'bold');
xlabel('Horn Length cut [m]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Horn Radius cut [m]', 'FontSize', 12, 'FontWeight', 'bold');
title('Optimized Antenna with respect to Return loss', 'FontSize', 12, 'FontWeight', 'bold');
grid on;