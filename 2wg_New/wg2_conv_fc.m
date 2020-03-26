A = load('SPP_conv.mat');

F = 4e9:0.5e9:21e9;

R = 0.0405319403216/2.1;

f =  fc(R, 1, 1);
    

for i = 23:30
    freq = f(i) + eps;
    j = find(abs(F - freq) < 1e8);
    if isempty(j)
        j = find(abs(F - freq) < 4e8);
        j = min(j);
    elseif isempty(j)
        j = find(abs(F - freq) < 1e9);
        j = min(j);
    end
    Spp(i) = A.Spp(i).Spp_i(j);
end

