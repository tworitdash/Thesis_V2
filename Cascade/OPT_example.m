

% [STT, STR, SRT, SRR, N] = GSM_N(R, L, er, mur, F, k);
c0 = 3e8;
F = 5e9;

lamb = c0./F;

Len = 5e-2;

fun = @(x) GSM_N_opt(x(1), x(2), x(3), Len, F, 20);

% [x, fval] = ga(fun, 2, [], [], [], [], [2e-2 4e-2], [4e-2 8e-2]);

x0 = [2e-2, 4e-2, round(Len./(lamb./10))];

% a = fun(x0)

A = [];
b = [];
Aeq = [];
beq = [];
lb = [2e-2, 4e-2, 2];
ub = [4e-2, 8e-2, round(Len./(lamb./40))];

tic;

[x, fval, ef, output, lambda] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub);

toc;



