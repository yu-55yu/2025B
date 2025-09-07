clc;
clear;
close all;

% --- �ļ��Ͳ�������
file1.path = '����1.xlsx';
file1.angle = 10; % ����� (��)
file2.path = '����2.xlsx';
file2.angle = 15; % ����� (��)

% --- ���ȫ������
config.waveN_fit_min = 1500; % �����ʼ���� (cm??)
config.n1_init = 2.58;       % ���Ӳ������ʳ�ֵ (����FFT������)
config.n2_real_init = 2.55;  % �ĵ������ʳ�ֵ
config.k1_A_init = 0.001;    % ���Ӳ�����ϵ������A��ֵ
config.k1_B_init = 2.0;      % ���Ӳ�����ϵ������B��ֵ

% Sellmeierģ�� B1, B2, B3, C1, C2, C3 ������ֵ
config.selParam_init = [5.5, 0.2, 0.05, 0.027, 100, 0.01]; 



% --- ���� 1: ʹ��FFT�����ȳ�ֵ
disp('=====================================================');
disp('���� 1: �����ȳ�ֵ');
data1 = readmatrix(file1.path);
thk_init1 = fft_thk_estimate(data1(:,1), data1(:,2)/100, config.n1_init, config.waveN_fit_min, file1.angle);

data2 = readmatrix(file2.path);
thk_init2 = fft_thk_estimate(data2(:,1), data2(:,2)/100, config.n1_init, config.waveN_fit_min, file2.angle);

% ʹ�������Ƕȹ�������ƽ��ֵ��Ϊȫ�ֳ�ֵ
config.thk_init = mean([thk_init1, thk_init2]);
fprintf('�ļ�1 (10��) FFT������: %.2f ��m\n', thk_init1);
fprintf('�ļ�2 (15��) FFT������: %.2f ��m\n', thk_init2);
fprintf('<strong>ƽ����ȳ�ֵ: %.2f ��m</strong>\n', config.thk_init);
disp('=====================================================');

% --- ���� 2: ��ÿ���ļ����ж�����ȫ�����
disp(' ');
disp('���� 2: ��ÿ���ļ��������');
disp('--> ��ʼ������1 (10��)...');
[result1] = process(data1, file1.angle, config, file1.path);
disp('--> ��ʼ������2 (15��)...');
[result2] = process(data2, file2.angle, config, file2.path);
disp('=====================================================');


% --- ���� 3: �����ͱȽ������Ƕȵ���Ͻ��
disp(' ');
disp('���� 3: �ۺϷ������');
analyze_res(result1, result2, file1.angle, file2.angle);
disp('=====================================================');




%% ���ݴ�������������
function [result] = process(data, angle, config, output_filepath)
epsilon= 11.7;
    waveNum = data(:, 1);       % ���� (cm??)
    R = data(:, 2) / 100;       % ������ (0-1)
    waveLen_full = 10000 ./ waveNum; % ���� (��m)

    % �����趨�Ĳ�����Χɸѡ������ϵ�����
    filter = waveNum >= config.waveN_fit_min;
    waveNum_fit = waveNum(filter);
    waveLen_fit = 10000 ./ waveNum_fit; 

    R_fit = R(filter);
    
    % ����˳��: [���, �ĵ�n2, Sellmeier(6��), k����(2��)]
    x0 = [config.thk_init, config.n2_real_init, config.selParam_init, config.k1_A_init, config.k1_B_init];
    lb = [config.thk_init*0.8, 2.0, 0.01, 0.001, 0.0001, 0.0001, 0.1, 0.001, 0, 0];
    ub = [config.thk_init*1.2, 3.5, 20, 10, 5, 2, 150, 20, 0.01, 4];

    model_type = 'real_substrate'; 
    [x_optimal, R_squared_fit] = global_fit(x0, lb, ub,  R_fit,waveLen_fit, angle, epsilon ,model_type);

    result.thk = x_optimal(1);
    result.n2 = x_optimal(2);
    result.selParam = x_optimal(3:8);
    result.k1Param = x_optimal(9:10);
    result.n1_complex_full = cal_n_sellmeier(waveLen_full, result.selParam, result.k1Param);

    fprintf('    ���պ��: %.2f ��m\n', result.thk);
    fprintf('    �������Ӳ� n1 (�� 6��m): %.3f + %.4fi\n', real(interp1(waveLen_full, result.n1_complex_full, 6)), imag(interp1(waveLen_full, result.n1_complex_full, 6)));
    fprintf('    ���ճĵ� n2: %.3f\n', result.n2);
    fprintf('    ��������������Ŷ� R?: %.4f\n', R_squared_fit);

    theta0_rad = angle * pi / 180;
    R_fit_full = compute_R(result.thk, result.n1_complex_full, result.n2, waveNum, theta0_rad,0);
    
    try
        R_fit_per = R_fit_full * 100;
        header = {'��Ϸ�����'};
        writecell(header, output_filepath, 'Sheet', 1, 'Range', 'C1');
        writematrix(R_fit_per, output_filepath, 'Sheet', 1, 'Range', 'C2');
        fprintf('    �ɹ�����Ͻ�����浽: %s\n', output_filepath);
    catch ME
        fprintf('    ����Excel�ļ�ʱ����: %s\n', ME.message);
    end

    plot_fit_res(waveNum, R, R_fit_full, angle, config.waveN_fit_min);
    plot_optical_constants(waveLen_full, waveNum, result.n1_complex_full, '���Ӳ� (Sellmeier)');
end


function analyze_res(result1, result2, angle1, angle2)
% ANALYZE_RES �����ͱȽ������Ƕȵ���Ͻ��
%
% ����:
%   result1 - ��һ���Ƕȵ���Ͻ���ṹ��
%   result2 - �ڶ����Ƕȵ���Ͻ���ṹ��
%   angle1  - ��һ���Ƕȵ�ֵ (��)
%   angle2  - �ڶ����Ƕȵ�ֵ (��)

% --- 1. �ı���������������бȽϹؼ����� ---

fprintf('<strong>--- �ؼ���������Ա� ---</strong>\n');

% �Ƚ����Ӳ��� (thk)
thk1 = result1.thk;
thk2 = result2.thk;
thk_diff_abs = abs(thk1 - thk2);
thk_diff_rel = thk_diff_abs / mean([thk1, thk2]) * 100;
fprintf('���Ӳ���:\n');
fprintf('  %d�� ��Ͻ��: %.3f ��m\n', angle1, thk1);
fprintf('  %d�� ��Ͻ��: %.3f ��m\n', angle2, thk2);
fprintf('  -> ���Բ���: %.3f ��m\n', thk_diff_abs);
fprintf('  -> ��Բ���: %.2f %%\n', thk_diff_rel);
if thk_diff_rel < 1.0
    fprintf('  ����: ��Ƚ�����и߶�һ���ԡ�\n\n');
else
    fprintf('  ����: ��Ƚ������һ�����죬���顣\n\n');
end

% �Ƚϳĵ������� (n2)
n2_1 = result1.n2;
n2_2 = result2.n2;
n2_diff_abs = abs(n2_1 - n2_2);
n2_diff_rel = n2_diff_abs / mean([n2_1, n2_2]) * 100;
fprintf('�ĵ������� (n2):\n');
fprintf('  %d�� ��Ͻ��: %.3f\n', angle1, n2_1);
fprintf('  %d�� ��Ͻ��: %.3f\n', angle2, n2_2);
fprintf('  -> ���Բ���: %.3f\n', n2_diff_abs);
fprintf('  -> ��Բ���: %.2f %%\n', n2_diff_rel);
if n2_diff_rel < 2.0
    fprintf('  ����: �ĵ������ʽ��һ�������á�\n\n');
else
    fprintf('  ����: �ĵ������ʽ������һ�����죬���顣\n\n');
end
end