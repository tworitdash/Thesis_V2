function [Eth, Eph] = Feed_FF_Superposition(ModeNumberAper, Gamma, Dm, theta, phi, F, er, mur, R, Transmission_sum)

Eth = zeros(size(theta));
Eph = zeros(size(theta));

for o = 1:ModeNumberAper
    
    HigherModes = o+1:1:o+20;
    
    [Eth_o, Eph_o] = FF_apertureFSCir(o, length(HigherModes)+1, [1, Dm(o, :)], Gamma(o), theta, phi, F, er(end), mur(end), R(end));
    Eth = Eth + Eth_o .* Transmission_sum(o);
    Eph = Eph + Eph_o .* Transmission_sum(o);

end