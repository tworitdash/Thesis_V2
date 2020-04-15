


S_14 = read(rfdata.data,'../../../feko/Cone_S.s48p'); 
s_14 = extract(S_14,'S_PARAMETERS');

str = zeros(38, 10);

for i = 1:38
    for j = 1:10
        str(i, j) = s_14(i, 38 + j);
    end
end

stt = zeros(38, 38);

for i = 1:38
    for j = 1:38
        stt(i, j) = s_14(i, j);
    end
end

% S = load('Str5_ratio_1_modes_20.mat');
% 
% 
% s_obs = squeeze(S.STR(21, 1:13, 1:10));

s_obs = squeeze(STR(1, :, :));

angle_diff = (angle(str) - angle(s_obs)) * 180/pi;

abs_diff = abs(str) - abs(s_obs);

real_diff = abs(real(str) - real(s_obs));
imag_diff = abs(imag(str) - imag(s_obs));


real_diff_abs = (abs(real(str)) - abs(real(s_obs)));
imag_diff_abs = (abs(imag(str)) - abs(imag(s_obs)));

figure;

surface(1:1:10, 1:1:38, (angle_diff), 'LineWidth', 2); shading flat;

grid on;

figure;

surface(1:1:10, 1:1:38, (abs_diff), 'LineWidth', 2); shading flat;

grid on;
% 
% figure;
% 
% surface(1:1:10, 1:1:13, (real_diff), 'LineWidth', 2); shading flat;
% 
% grid on;
% 
% figure;
% 
% surface(1:1:10, 1:1:13, (imag_diff), 'LineWidth', 2); shading flat;
% 
% grid on;
% 
% figure;
% 
% surface(1:1:10, 1:1:13, (real_diff_abs), 'LineWidth', 2); shading flat;
% 
% grid on;
% 
% figure;
% 
% surface(1:1:10, 1:1:13, (imag_diff_abs), 'LineWidth', 2); shading flat;
% 
% grid on;
%  
% 
% figure;
% 
% surface(1:1:10, 1:1:13, angle(str)*180/pi, 'LineWidth', 2); shading flat;
% 
% grid on;
% 
% figure;
% 
% surface(1:1:10, 1:1:13, angle(s_obs)*180/pi, 'LineWidth', 2); shading flat;
% 
% grid on;

figure;

surface(1:1:10, 1:1:38, abs(str), 'LineWidth', 2); shading flat;

grid on;

figure;

surface(1:1:20, 1:1:78, db(abs(s_obs)), 'LineWidth', 2); shading flat;

grid on;