clc;
clear;
close all;

file1.path = '����1.xlsx';
file1.angle = 10; % ����� (��)
file2.path = '����2.xlsx';
file2.angle = 15; % ����� (��)

config.waveN_fit_min = 1500; % �����ʼ���� (cm-1)
config.n1_init = 2.58;       % ���Ӳ������ʳ�ֵ (����FFT������)
config.n2_real_init = 2.55;  % �ĵ������ʳ�ֵ
config.k1_A_init = 0.001;    % ���Ӳ�����ϵ������A��ֵ
config.k1_B_init = 2.0;      % ���Ӳ�����ϵ������B��ֵ

% selģ�� B1, B2, B3, C1, C2, C3 ������ֵ
config.selParam_init = [5.5, 0.2, 0.05, 0.027, 100, 0.01];



disp('FFT�����ȳ�ֵ');
data1 = readmatrix(file1.path);
thk_init1 = fft_thk_estimate(data1(:,1), data1(:,2)/100, config.n1_init, config.waveN_fit_min, file1.angle);
data2 = readmatrix(file2.path);
thk_init2 = fft_thk_estimate(data2(:,1), data2(:,2)/100, config.n1_init, config.waveN_fit_min, file2.angle);
config.thk_init = mean([thk_init1, thk_init2]);
fprintf('�ļ�1 (10��) FFT������: %.2f ��m\n', thk_init1);
fprintf('�ļ�2 (15��) FFT������: %.2f ��m\n', thk_init2);
fprintf('<strong>ƽ����ȳ�ֵ: %.2f ��m</strong>\n', config.thk_init);


disp(' ');
disp('10��');
[res1] = process(data1, file1.angle, config, file1.path);
disp('15��');
[res2] = process(data2, file2.angle, config, file2.path);


disp(' ');
disp('������');
analyze_res(res1, res2, file1.angle, file2.angle);




%% ���ݴ�������������
function [res] = process(data, angle, config, output_filepath)
epsilon= 11.7;
waveNum = data(:, 1);       % ���� (cm-1)
R = data(:, 2) / 100;       % ������ (0-1)
waveLen = 10000 ./ waveNum; % ���� (��m)

filter = waveNum >= config.waveN_fit_min;
waveNum_fit = waveNum(filter);
waveLen_fit = 10000 ./ waveNum_fit;

R_fit = R(filter);

% ���, �ĵ�n2, sel(6), k����(2)
x0 = [config.thk_init, config.n2_real_init, config.selParam_init, config.k1_A_init, config.k1_B_init];
lb = [config.thk_init*0.8, 2.0, 0.01, 0.001, 0.0001, 0.0001, 0.1, 0.001, 0, 0];
ub = [config.thk_init*1.2, 3.5, 20, 10, 5, 2, 150, 20, 0.01, 4];

model_type = 'real_substrate';
[x_optimal, R_squared_fit] = global_fit(x0, lb, ub,  R_fit,waveLen_fit, angle, epsilon ,model_type);

res.thk = x_optimal(1);
res.n2 = x_optimal(2);
res.selParam = x_optimal(3:8);
res.k1Param = x_optimal(9:10);
res.n1_complex = cal_n_sellmeier(waveLen, res.selParam, res.k1Param);

fprintf('    ���պ��: %.2f ��m\n', res.thk);
fprintf('    �������Ӳ� n1 (�� 6��m): %.3f + %.4fi\n', real(interp1(waveLen, res.n1_complex, 6)), imag(interp1(waveLen, res.n1_complex, 6)));
fprintf('    ���ճĵ� n2: %.3f\n', res.n2);
fprintf('    ��������������Ŷ� R?: %.4f\n', R_squared_fit);

theta0_rad = angle * pi / 180;


R_fit_full = compute_R(res.thk, res.n1_complex, res.n2, waveNum, theta0_rad,2);


R_fit_per = R_fit_full * 100;
header = {'��Ϸ�����'};
writecell(header, output_filepath, 'Sheet', 1, 'Range', 'C1');
writematrix(R_fit_per, output_filepath, 'Sheet', 1, 'Range', 'C2');


plot_fit_res(waveNum, R, R_fit_full, angle, config.waveN_fit_min);
plot_optical_constants(waveLen, waveNum, res.n1_complex, '���Ӳ� (sel)');
end


function analyze_res(res1, res2, angle1, angle2)
% ANALYZE_RES �����ͱȽ������Ƕȵ���Ͻ��
%
% ����:
%   res1 - ��һ���Ƕȵ���Ͻ���ṹ��
%   res2 - �ڶ����Ƕȵ���Ͻ���ṹ��
%   angle1  - ��һ���Ƕȵ�ֵ (��)
%   angle2  - �ڶ����Ƕȵ�ֵ (��)


fprintf('<strong> �ؼ���������Ա� </strong>\n');

% �Ƚ����Ӳ��� (thk)
thk1 = res1.thk;
thk2 = res2.thk;
thk_diff_abs = abs(thk1 - thk2);
thk_diff_rel = thk_diff_abs / mean([thk1, thk2]) * 100;
fprintf('���Ӳ���:\n');
fprintf('  %d�� ��Ͻ��: %.3f ��m\n', angle1, thk1);
fprintf('  %d�� ��Ͻ��: %.3f ��m\n', angle2, thk2);
fprintf('  -> ���Բ���: %.3f ��m\n', thk_diff_abs);
fprintf('  -> ��Բ���: %.2f %%\n', thk_diff_rel);
if thk_diff_rel < 10.0
    fprintf('  ����: ��Ƚ�����и߶�һ���ԡ�\n\n');
else
    fprintf('  ����: ��Ƚ������һ�����죬���顣\n\n');
end

% �Ƚϳĵ������� (n2)
n2_1 = res1.n2;
n2_2 = res2.n2;
n2_diff_abs = abs(n2_1 - n2_2);
n2_diff_rel = n2_diff_abs / mean([n2_1, n2_2]) * 100;
fprintf('�ĵ������� (n2):\n');
fprintf('  %d�� ��Ͻ��: %.3f\n', angle1, n2_1);
fprintf('  %d�� ��Ͻ��: %.3f\n', angle2, n2_2);
fprintf('  -> ���Բ���: %.3f\n', n2_diff_abs);
fprintf('  -> ��Բ���: %.2f %%\n', n2_diff_rel);
if n2_diff_rel < 10.0
    fprintf('  ����: �ĵ������ʽ��һ�������á�\n\n');
else
    fprintf('  ����: �ĵ������ʽ������һ�����죬���顣\n\n');
end
end