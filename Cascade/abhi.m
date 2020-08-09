a = 2201;
b = 5800;

x = a:300:b;

r = a + round((b - a)*rand(304, 1));

for i = 2:length(x)
    j = x(i);
    k = x(i - 1);
    
    f = find((r > k) & (r < j));
    Xn(i - 1).n = r(f);
    
    num_in_the_interval(1, i - 1) = size(f, 1);
end

max_col_size = max(num_in_the_interval);

out = zeros(length(x) - 1, max_col_size);

for m = 1:size(out, 1)-1
    out(m, 1:length(Xn(m).n)) = (Xn(m).n).';
end

