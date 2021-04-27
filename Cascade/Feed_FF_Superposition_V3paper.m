function [Eth, Eph, Eco, Exp, CO, XP] = Feed_FF_Superposition_V3paper(ModeNumberAper, Gamma, Dm, theta, phi, F, er, mur, R, Transmission_sum, HM)

Eth = zeros(size(theta));
Eph = zeros(size(theta));
Eco = zeros(size(theta));
Exp = zeros(size(theta));
CO = zeros(size(theta));
XP = zeros(size(theta));


for o = 1:ModeNumberAper
    
    HigherModes = ModeNumberAper+1:1:HM;
    
    [Eth_o, Eph_o, Eco_o, Exp_o, CO_o, XP_o] = FF_apertureFSCir4_paper([o, HigherModes], [1, Dm(o, ModeNumberAper+1:end)], Gamma(o), theta, phi, F, er(end), mur(end), R(end), R(1));
    Eth = Eth + Eth_o .* Transmission_sum(o);
    Eph = Eph + Eph_o .* Transmission_sum(o);
    Eco = Eco + Eco_o .* Transmission_sum(o);
    Exp = Exp + Exp_o .* Transmission_sum(o);
    CO = CO + CO_o .* Transmission_sum(o);
    XP = XP + XP_o .* Transmission_sum(o);
end