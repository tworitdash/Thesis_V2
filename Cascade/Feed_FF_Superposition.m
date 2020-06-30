function [Eth, Eph, Eco, Exp, CO, XP] = Feed_FF_Superposition(ModeNumberAper, Gamma, Dm, theta, phi, F, er, mur, R, Transmission_sum)

Eth = zeros(size(theta));
Eph = zeros(size(theta));
Eco = zeros(size(theta));
Exp = zeros(size(theta));
CO = zeros(size(theta));
XP = zeros(size(theta));


for o = 1:ModeNumberAper
    
    HigherModes = o+1:1:20;
    
    [Eth_o, Eph_o, Eco_o, Exp_o, CO, XP] = FF_apertureFSCir(o, length(HigherModes)+1, [1, Dm(o, :)], Gamma(o), theta, phi, F, er(end), mur(end), R(end));
    Eth = Eth + Eth_o .* Transmission_sum(o);
    Eph = Eph + Eph_o .* Transmission_sum(o);
    Eco = Eco + Eco_o .* Transmission_sum(o);
    Exp = Exp + Exp_o .* Transmission_sum(o);
    CO = CO + CO .* Transmission_sum(o);
    XP = XP + XP .* Transmission_sum(o);
end