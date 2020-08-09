function I2rho = intphicos(a, b, pm, rm)

if pm == rm
    I2rho = (1./2) .* ( (b - a) + (1./(2.*pm)) .* (sin(2.*pm.*b) - sin(2.*pm.*a)));
else

I2rho = (1./2) .* ( (1./(pm - rm) .* ( sin(b.*(pm - rm)) - sin(a .*(pm - rm)))) + ...
(1./(pm + rm) .* ( sin(b.*(pm + rm)) - sin(a .*(pm + rm)))) );
end

end