
N = 25;

for i = 1:N
    for j = 1:i
            Sper_1(sum(0:i-1) + j).mode = "TE";
            Sper_1(sum(0:i-1) + j).pol = "even";
            Sper_1(sum(0:i-1) + j).n = i;
            Sper_1(sum(0:i-1) + j).m = j;
    end
end

for i = 1:N
    for j = 1:i
            Sper_2(sum(0:i-1) + j).mode = "TE";
            Sper_2(sum(0:i-1) + j).pol = "odd";
            Sper_2(sum(0:i-1) + j).n = i;
            Sper_2(sum(0:i-1) + j).m = j;
    end
end

for i = 1:N
    for j = 1:i
            Sper_3(sum(0:i-1) + j).mode = "TM";
            Sper_3(sum(0:i-1) + j).pol = "even";
            Sper_3(sum(0:i-1) + j).n = i;
            Sper_3(sum(0:i-1) + j).m = j;
    end
end

for i = 1:N
    for j = 1:i
            Sper_4(sum(0:i-1) + j).mode = "TM";
            Sper_4(sum(0:i-1) + j).pol = "odd";
            Sper_4(sum(0:i-1) + j).n = i;
            Sper_4(sum(0:i-1) + j).m = j;
    end
end


Sper = [Sper_1, Sper_2, Sper_3, Sper_4];

save('Sper.mat', 'Sper');