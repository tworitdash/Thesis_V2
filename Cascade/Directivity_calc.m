c0 = 3e8;
N = 100;
F = 5e9;

lamb = c0/F;

Len = 6 .* lamb;

R1 = 2e-2;
Rend = 2 .* lamb;

R_test = linspace(R1, Rend, N);


problem2.objective = @(x) MinXP_Goal_V2L_fields([x(1:N)], x(N+1), F, 0);

IP = [R_test Len];

lb = zeros(N+1, 1);
ub = zeros(N+1, 1);

lb(1) = R1;
ub(1) = R1 + R1/5;

for i =  2:N
    lb(i) = ub(i - 1) - R_test(1)/5;
    ub(i) = R_test(i) + R_test(1)/5;
end
ub(N + 1) = 8 .* lamb;
lb(N + 1) = 5 .* lamb; 

zeta = 120 * pi;

[Eth_Ip, Eph_Ip, Eco_Ip, Exp_Ip, CO, XP, the, ph] = problem2.objective(IP);
% Eco_Ip = Eth_Ip(46, :) .* cos(the(1, :)) + Eph_Ip(46, :);
% Exp_Ip = Eth_Ip(46, :) .* cos(the(1, :)) - Eph_Ip(46, :);

E_Ip = sqrt(abs(Eth_Ip).^2 + abs(Eph_Ip).^2);
U_feedIp = abs(E_Ip).^2./(2 .* zeta);
P_radIp_i = U_feedIp(:, 91:end) .* sin(the(1, 91:end)) .* pi./180 .* pi./180;
P_radIp = sum(sum(P_radIp_i));
D_Ip = 4 .* pi .* U_feedIp ./ P_radIp;

% [Eco_Lb, Exp_Lb, ~, ~] = problem2.objective(lb);
% E_Lb = sqrt(abs(Eco_Lb).^2 + abs(Exp_Lb).^2);
% U_feedLb = abs(E_Lb).^2./(2 .* zeta);
% P_radLb_i = U_feedLb(:, 91:end) .* sin(the(1, 91:end)) .* pi./180 .* pi./180;
% P_radLb = sum(sum(P_radLb_i));
% D_Lb = 4 .* pi .* U_feedLb ./ P_radLb;
% 
% [Eco_Ub, Exp_Ub, ~, ~] = problem2.objective(ub);
% E_Ub = sqrt(abs(Eco_Ub).^2 + abs(Exp_Ub).^2);
% U_feedUb = abs(E_Ub).^2./(2 .* zeta);
% P_radUb_i = U_feedUb(:, 91:end) .* sin(the(1, 91:end)) .* pi./180 .* pi./180;
% P_radUb = sum(sum(P_radUb_i));
% D_Ub = 4 .* pi .* U_feedUb ./ P_radUb;
%%

Exp_Ip = 1./sqrt(2) .* (Eth_Ip(46, :) - Eph_Ip(46, :));
Eco_Ip = 1./sqrt(2) .* (Eth_Ip(46, :) + Eph_Ip(46, :));

figure(5);
hold on;
plot(the(1, :)*180/pi, db(abs(Eco_Ip)./max(abs(Eco_Ip))), 'LineWidth', 2);
hold on;
plot(the(1, :)*180/pi, db(abs(Exp_Ip)./max(abs(Eco_Ip))), 'LineWidth', 2);
grid on;
ylim([-50 0]);

hold on;
plot(the(1, :)*180/pi, db(abs(CO(1, :))./max(abs(CO(1, :)))), 'LineWidth', 2);
hold on;
plot(the(1, :)*180/pi, db(abs(XP(46, :))./max(abs(CO(1, :)))), 'LineWidth', 2);
grid on;
ylim([-50 0]);


figure(10);
hold on;
plot(the(1, :)*180/pi, db(abs(Eth_Ip(46, :))./max(abs(Eth_Ip(46, :)))), 'LineWidth', 2);
hold on;
plot(the(1, :)*180/pi, db(abs(Eph_Ip(46, :))./max(abs(Eph_Ip(46, :)))), 'LineWidth', 2);
grid on;
ylim([-50 0]);

figure;
plot(the(1, :)*180/pi, db(abs(Eco_Ip(1, :))./max(abs(Eco_Ip(1, :)))), 'b', 'LineWidth', 2);
hold on;
plot(the(1, :)*180/pi, db(abs(Exp_Ip(46, :))./max(abs(Eco_Ip(1, :)))), 'LineWidth', 2, 'color', [0, 1, 0.2]);
hold on;
ylim([-50 0]);
% plot(the(1, :)*180/pi, db(abs(Eco_Lb(1, :))./max(abs(Eco_Lb(1, :)))), 'r', 'LineWidth', 2);
% hold on;
% plot(the(1, :)*180/pi, db(abs(Exp_Lb(46, :))./max(abs(Eco_Lb(1, :)))), 'LineWidth', 2, 'color', [1, 0, 0.2]);
% hold on;
% plot(the(1, :)*180/pi, db(abs(Eco_Ub(1, :))./max(abs(Eco_Ub(1, :)))), 'g', 'LineWidth', 2);
% hold on;
% plot(the(1, :)*180/pi, db(abs(Exp_Ub(46, :))./max(abs(Eco_Ub(1, :)))), 'LineWidth', 2, 'color', [0.2, 0, 1]);
% grid on;
% hold on;
% plot(the(1, :)*180/pi, db(abs(Eco_Ub(1, :))./max(abs(Eco_Ub(1, :)))), 'g', 'LineWidth', 2);

figure;
hold on;
plot(the(1, :)*180/pi, db(abs(E_Ip(1, :))./max(abs(E_Ip(1, :)))), 'LineWidth', 2);
hold on;
plot(the(1, :)*180/pi, db(abs(E_Ip(91, :))./max(abs(E_Ip(91, :)))), 'LineWidth', 2);
grid on;
ylim([-50 0]);

% figure;
% plot(the(1, :)*180/pi, db(D_Ip(1, :))/2, 'LineWidth', 2);
% hold on;
% plot(the(1, :)*180/pi, db(D_Lb(1, :))/2, 'LineWidth', 2);
% hold on;
% plot(the(1, :)*180/pi, db(D_Ub(1, :))/2, 'LineWidth', 2);
% grid on;


%% 

l1 = lamb/4;
lIP = [l1 ones(1, N - 1).* Len/N];
llb = [l1 ones(1, N - 1).* lb(end)/N];
lub = [l1 ones(1, N - 1).* ub(end)/N];
for i = 1:N
    lIP_axis(i) = sum(lIP(1:i));
    llb_axis(i) = sum(llb(1:i));
    lub_axis(i) = sum(lub(1:i));
end

figure;
plot([lIP_axis], R_test, 'b', 'LineWidth', 2); hold on; plot([lIP_axis], -R_test, 'b', 'LineWidth', 2);
hold on;
plot([llb_axis], lb(1:N), 'r', 'LineWidth', 2); hold on; plot([llb_axis], -lb(1:N), 'r', 'LineWidth', 2);
hold on;
plot([lub_axis], ub(1:N), 'g', 'LineWidth', 2); hold on; plot([lub_axis], -ub(1:N), 'g', 'LineWidth', 2);

grid on;

figure(6);
hold on;
plot(the(1, :)*180/pi, db(abs(Eth_Ip(46, :))), 'LineWidth', 2);
hold on;
plot(the(1, :)*180/pi, db(abs(Eph_Ip(46, :))), 'LineWidth', 2);
grid on;
ylim([-50 0]);
