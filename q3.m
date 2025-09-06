clc;
clear;
close all;

file1.path = '����3.xlsx';
file1.angle = 10;
file2.path = '����4.xlsx';
file2.angle = 15;

% ��ϵĲ�����Χ
config.waveN_fit_min = 400;

%% �����趨
config.n1_init = 3.48; % ���Ӳ������ʳ�ֵ

% ���Ӳ����ģ�� n(lambda) = A + B / lambda^2 ����ϳ�ֵ
config.cauchyParam_init = [3.42, 0.05]; 

% ���Ӳ�����ģ�� k1 = A * lambda^B
config.k1_A_init = 1e-5; 
config.k1_B_init = 2.0;

% �ĵ�Drudeģ�Ͳ�����ֵ [nu_p (cm^-1), Gamma (cm^-1)]
% nu_p ~ ��������Ƶ�ʣ������Ũ���й�; Gamma ~ ����ϵ��
config.drudeParam_init = [2000, 200]; 
config.epsilon_inf = 11.7; % ��ĸ�Ƶ��糣�� (�̶�ֵ)



disp('�����ȳ�ֵ');
data1 = readmatrix(file1.path);
thk_init1 = fft_thk_estimate(data1(:,1), data1(:,2)/100, config.n1_init, 2000, file1.angle); % FFT�����ڸ�Ƶ������
data2 = readmatrix(file2.path);
thk_init2 = fft_thk_estimate(data2(:,1), data2(:,2)/100, config.n1_init, 2000, file2.angle);
config.thk_init = mean([thk_init1, thk_init2]);
fprintf('�ļ�3 (10��) FFT������: %.2f ��m\n', thk_init1);
fprintf('�ļ�4 (15��) FFT������: %.2f ��m\n', thk_init2);
fprintf('ƽ����ȳ�ֵ: %.2f ��m\n', config.thk_init);
fprintf('\n');

disp('--- ��ʼ������3 (10��) ---')
[result1] = process_file(data1, file1.angle, config);
fprintf('\n');
disp('--- ��ʼ������4 (15��) ---')
[result2] = process_file(data2, file2.angle, config);

analyze_results(result1, result2, file1.angle, file2.angle);

%% ���Ĵ�����
function [result] = process_file(data, incident_angle, config)
    waveNum = data(:, 1);
    R = data(:, 2) / 100;
    waveLen_full = 10000 ./ waveNum;
    filter = waveNum > config.waveN_fit_min; 
    waveNum_fit = waveNum(filter);
    R_fit = R(filter);
    waveLen_fit = waveLen_full(filter);
    
    %% ��������
    % ���(1), ����(2), ���Ӳ�����(2), �ĵ�Drude(2)
    x0 = [config.thk_init, config.cauchyParam_init, config.k1_A_init, config.k1_B_init, config.drudeParam_init];
    
    % thk,    cauchy A, B,    k1_A, B,   drude nu_p, Gamma
    lb = [config.thk_init*0.9, 3.4, 0,    0,    0,   100,  10];
    ub = [config.thk_init*1.1, 3.5, 0.2,  1e-3, 4,   4000, 1000];

    [x_optimal, R_squared_fit] = global_fit_all_parameters(x0, lb, ub, waveNum_fit, R_fit, waveLen_fit, incident_angle, config.epsilon_inf);

    %% ���Ż�����зֽ����
    result.thk = x_optimal(1);
    result.cauchyParam = x_optimal(2:3);
    result.k1Param = x_optimal(4:5);
    result.drudeParam = x_optimal(6:7);
    
    % ����ȫ�׶εĹ�ѧ����
    result.n1_complex_full = calculate_n1_complex(waveLen_full, result.cauchyParam, result.k1Param);
    result.n2_complex_full = calculate_n2_complex_drude(waveNum, result.drudeParam, config.epsilon_inf);

    fprintf('  ���պ��: %.2f ��m\n', result.thk);
    fprintf('  �������Ӳ� ?0?91 (6��m): %.3f + %.4fi\n', real(interp1(waveLen_full, result.n1_complex_full, 6)), imag(interp1(waveLen_full, result.n1_complex_full, 6)));
    fprintf('  ���ճĵ�Drude���� ��p: %.1f cm?6?3?0?1, ��: %.1f cm?6?3?0?1\n', result.drudeParam(1), result.drudeParam(2));
    fprintf('  ��������������Ŷ� R?0?5: %.4f\n', R_squared_fit);

    plot_fitting_results(waveNum, R, result.n1_complex_full, result.n2_complex_full, result.thk, incident_angle, config);
    plot_refractive_index(waveLen_full, waveNum, result.n1_complex_full, '���Ӳ�');
    plot_refractive_index(waveLen_full, waveNum, result.n2_complex_full, '�ĵ� (Drudeģ��)');
end

%% ȫ���Ż�����
function [x_optimal, R_squared] = global_fit_all_parameters(x0, lb, ub, waveNum_fit, R_fit, waveLen_fit, incident_angle, epsilon_inf)
    theta0_rad = incident_angle * pi / 180;
    % ��epsilon_inf����ģ��
    model_func = @(x, k) model_R(x, k, waveLen_fit, theta0_rad, epsilon_inf);
    options = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxIterations', 1000, ...
        'FunctionTolerance', 1e-9, 'StepTolerance', 1e-10, ...
        'Algorithm', 'trust-region-reflective', 'UseParallel', true);
    [x_optimal, ~] = lsqcurvefit(model_func, x0, waveNum_fit, R_fit, lb, ub, options);
    R_fitted = model_func(x_optimal, waveNum_fit);
    SS_tot = sum((R_fit - mean(R_fit)).^2);
    SS_res = sum((R_fit - R_fitted).^2);
    R_squared = 1 - SS_res / SS_tot;
end

%% FFT ��ȹ���
function thk = fft_thk_estimate(waveNum, R, n_avg, waveN_fit_min, theta0_deg)
    if nargin < 5
        theta0_deg = 0; 
    end 
    theta0_rad = theta0_deg * pi / 180; 
    cos_theta1 = real(sqrt(1 - (sin(theta0_rad) / n_avg)^2)); 
    filter = waveNum > waveN_fit_min; 
    waveNum_fft = waveNum(filter); 
    R_fft = R(filter); 
    R_ac = R_fft - mean(R_fft); 
    N = 2^nextpow2(8*length(waveNum_fft)); 
    k_uniform = linspace(min(waveNum_fft), max(waveNum_fft), N); 
    r_uniform = interp1(waveNum_fft, R_ac, k_uniform, 'pchip', 'extrap'); 
    window = hann(N); 
    r_windowed = r_uniform(:) .* window(:); fft_result = fft(r_windowed); 
    fft_power = abs(fft_result(1:N/2)).^2; dk = mean(diff(k_uniform)); 
    thk_axis = (0:N/2-1) * 10000 / (2 * n_avg * cos_theta1 * N * dk); 
    search_range = find(thk_axis > 5 & thk_axis < 200); 
    if isempty(search_range)
        thk = 20; 
        return; 
    end
    [~, max_idx_in_range] = max(fft_power(search_range)); 
    max_idx_rough = search_range(max_idx_in_range); 
    correctNum = 3; DatePower1 = 0; DatePower2 = 0; 
    for i = -correctNum:correctNum, idx = max_idx_rough + i; 
        if idx >= 1 && idx <= length(fft_power)
            power = fft_power(idx); 
            DatePower1 = DatePower1 + idx * power; 
            DatePower2 = DatePower2 + power; 
        end
    end

    if DatePower2 > 0
        f_corrected = DatePower1 / DatePower2; 
    else
        f_corrected = max_idx_rough;
    end

    thk = (f_corrected - 1) * 10000 / (2 * n_avg * cos_theta1 * N * dk);
end

%% �����������
function analyze_results(res1, res2, angle1, angle2)
    disp('--- ���ս���Ա� ---');
    fprintf('�ļ�3 (%d��) -> ���: %.2f ��m, ��p: %.1f, ��: %.1f\n', angle1, res1.thk, res1.drudeParam(1), res1.drudeParam(2));
    fprintf('�ļ�4 (%d��) -> ���: %.2f ��m, ��p: %.1f, ��: %.1f\n', angle2, res2.thk, res2.drudeParam(1), res2.drudeParam(2));
    thk1 = res1.thk; thk2 = res2.thk;
    thk_diff_percent = abs(thk1-thk2)/mean([thk1,thk2])*100;
    fprintf('�����ǶȲ�õĺ�ȷֱ�Ϊ %.2f ��m �� %.2f ��m����Բ���Ϊ %.2f%%��\n', thk1, thk2, thk_diff_percent);
    if thk_diff_percent < 5, fprintf('  -> ����: ���һ�������ã�����ɿ���\n'); fprintf('  -> �Ƽ����ֵ: %.2f �� %.2f ��m\n', mean([thk1, thk2]), std([thk1, thk2])); else, fprintf('  -> ����: ��Ȳ���ϴ�\n'); end
end

%% ����ģ�����ͼ����

function R = model_R(params, ~, waveLen, theta0, epsilon_inf) 
    thk = params(1);
    cauchyParam = params(2:3);
    k1Param = params(4:5);
    drudeParam = params(6:7);
    
    waveNum_model = 10000 ./ waveLen;
    n1_complex = calculate_n1_complex(waveLen, cauchyParam, k1Param);
    n2_complex = calculate_n2_complex_drude(waveNum_model, drudeParam, epsilon_inf);

    R = model_R_vectorized(thk, n1_complex, n2_complex, waveNum_model, theta0);
end

function n1_complex = calculate_n1_complex(waveLen, cauchyParam, k1Param)
    A = cauchyParam(1); B = cauchyParam(2);
    n_real = A + B ./ (waveLen.^2);
    k1_A = k1Param(1); k1_B = k1Param(2);
    k_imag = k1_A * waveLen.^k1_B;
    n1_complex = n_real + 1i * k_imag;
end

%% Drudeģ�ͺ���������ĵ׸������� 
function n2_complex = calculate_n2_complex_drude(waveNum, drudeParam, epsilon_inf)
    nu_p = drudeParam(1);   % ��������Ƶ�� (cm^-1)
    Gamma = drudeParam(2); % ����ϵ�� (cm^-1)
    
    nu = waveNum(:);
    
    epsilon_complex = epsilon_inf - (nu_p^2) ./ (nu.^2 + 1i * Gamma * nu);
    
    n2_complex = sqrt(epsilon_complex);
end

%% 
function R = model_R_vectorized(thk, n1_complex, n2_complex, waveNum, theta0)
    n0 = 1.0;
    n1 = n1_complex(:);
    n2 = n2_complex(:); 
    waveNum = waveNum(:);
    sin_theta1 = n0 * sin(theta0) ./ n1;
    cos_theta1 = sqrt(1 - sin_theta1.^2);
    sin_theta2 = n0 * sin(theta0) ./ n2;
    cos_theta2 = sqrt(1 - sin_theta2.^2); % n2�Ǹ���, cos_theta2Ҳ���Ǹ���
    
    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    r12_s = (n1.*cos_theta1 - n2.*cos_theta2) ./ (n1.*cos_theta1 + n2.*cos_theta2);
    r01_p = (n1*cos(theta0) - n0*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    r12_p = (n2.*cos_theta1 - n1.*cos_theta2) ./ (n2.*cos_theta1 + n1.*cos_theta2);
    
    delta = 4 * pi * n1 .* thk .* cos_theta1 .* waveNum / 10000;
    exp_term = exp(1i*delta);
    

    % �����
    % r_s_total = (r01_s + r12_s .* exp_term) ./ (1 + r01_s .* r12_s .* exp_term);
    % r_p_total = (r01_p + r12_p .* exp_term) ./ (1 + r01_p .* r12_p .* exp_term);

    % ˫��������
    r_s_total = r01_s + r12_s .* exp_term;
    r_p_total = r01_p + r12_p .* exp_term;

    R_s = abs(r_s_total).^2;
    R_p = abs(r_p_total).^2;
    R = (R_s + R_p) / 2;
    R = real(R(:));
end

function plot_refractive_index(waveLen, waveNum, n_complex_full, layer_name)
    figure('Name', ['��ϵõ���' layer_name '��������ɫɢ����'], 'Position', [100, 100, 1200, 500]);
    n_real = real(n_complex_full);
    k_imag = imag(n_complex_full);
    subplot(2,2,1); plot(waveLen, n_real, 'b-', 'LineWidth', 2); xlabel('���� (��m)'); ylabel(['������ʵ�� n (' layer_name ')']); title(['n vs ����']); grid on;
    subplot(2,2,2); plot(waveLen, k_imag, 'r-', 'LineWidth', 2); xlabel('���� (��m)'); ylabel(['����ϵ�� k (' layer_name ')']); title(['k vs ����']); grid on;
    subplot(2,2,3); plot(waveNum, n_real, 'b-', 'LineWidth', 2); xlabel('���� (cm^{-1})'); ylabel(['������ʵ�� n (' layer_name ')']); title(['n vs ����']); grid on; xlim([min(waveNum), max(waveNum)]);
    subplot(2,2,4); plot(waveNum, k_imag, 'r-', 'LineWidth', 2); xlabel('���� (cm^{-1})'); ylabel(['����ϵ�� k (' layer_name ')']); title(['k vs ����']); grid on; xlim([min(waveNum), max(waveNum)]);
end

function plot_fitting_results(waveNum, R, n1_complex_full, n2_complex_full, thk, incident_angle, config)
    theta0 = incident_angle * pi / 180;
    R_fitted = model_R_vectorized(thk, n1_complex_full, n2_complex_full, waveNum, theta0);
    R_squared_full = 1 - sum((R - R_fitted).^2) / sum((R - mean(R)).^2);
    figure('Name', ['��Ͻ������ (����� ' num2str(incident_angle) '��)'], 'Position', [150, 150, 1200, 600]);
    plot(waveNum, R*100, 'b.', 'MarkerSize', 5, 'DisplayName', 'ʵ������');
    hold on;
    plot(waveNum, R_fitted*100, 'r-', 'LineWidth', 2, 'DisplayName', '���ģ�� (Drude)');
    xlabel('���� (cm^{-1})'); ylabel('������ (%)');
    title_str = sprintf('��������Ͻ�� (��=%d��, d=%.2f ��m, ��_p=%.1f, ��=%.1f, ȫ��R?0?5=%.4f)', ...
        incident_angle, thk, config.drudeParam_init(1), config.drudeParam_init(2), R_squared_full);
    title(title_str);
    legend('Location', 'best'); grid on; xlim([min(waveNum), max(waveNum)]);
end