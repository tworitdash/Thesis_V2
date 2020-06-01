function [Dk] = Dkx(kx, k0, b)

dky = kx(2) - kx(1);

[kx, ky] = meshgrid(kx, kx);

kz = -1j .* sqrt(-(k0.^2 - kx.^2 - ky.^2));

Dk_int = 1./(kz) .* (sinc(ky .* b./2/pi)).^2;

Dk = sum(Dk_int, 1) .* dky;

end