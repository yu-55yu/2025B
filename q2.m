clc;
clear;
close all;

% --- �ļ��Ͳ�������
file1.path = '����1.xlsx';
file1.angle = 10; % ����� (��)
file2.path = '����2.xlsx';
file2.angle = 15; % ����� (��)

% --- ���ȫ������
config.waveN_fit_min = 1200; % �����ʼ���� (cm??)
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
    waveNum = data(:, 1);       % ���� (cm??)
    R = data(:, 2) / 100;       % ������ (0-1)
    waveLen_full = 10000 ./ waveNum; % ���� (��m)

    % �����趨�Ĳ�����Χɸѡ������ϵ�����
    filter = waveNum >= config.waveN_fit_min;
    waveNum_fit = waveNum(filter);
    R_fit = R(filter);
    
    % ����˳��: [���, �ĵ�n2, Sellmeier(6��), k����(2��)]
    x0 = [config.thk_init, config.n2_real_init, config.selParam_init, config.k1_A_init, config.k1_B_init];
    lb = [config.thk_init*0.8, 2.0, 0.1, 0.001, 0.0001, 0.0001, 0.1, 0.001, 0, 0];
    ub = [config.thk_init*1.2, 3.5, 20, 10, 5, 2, 150, 20, 0.01, 4];

    model_type = 'real_substrate'; 
    [x_optimal, R_squared_fit] = global_fit(x0, lb, ub, waveNum_fit, R_fit, angle, model_type);

    result.thk = x_optimal(1);
    result.n2 = x_optimal(2);
    result.selParam = x_optimal(3:8);
    result.k1Param = x_optimal(9:10);
    result.n1_complex_full = cal_n_sellmeier(waveLen_full, result.selParam, result.k1Param); % <--- ������������

    fprintf('    ���պ��: %.2f ��m\n', result.thk);
    fprintf('    �������Ӳ� n1 (�� 6��m): %.3f + %.4fi\n', real(interp1(waveLen_full, result.n1_complex_full, 6)), imag(interp1(waveLen_full, result.n1_complex_full, 6)));
    fprintf('    ���ճĵ� n2: %.3f\n', result.n2);
    fprintf('    ��������������Ŷ� R?: %.4f\n', R_squared_fit);

    theta0_rad = angle * pi / 180;
    R_fit_full = compute_R(result.thk, result.n1_complex_full, result.n2, waveNum, theta0_rad, 1);
    
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


