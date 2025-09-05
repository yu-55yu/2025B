clc;
clear;
close all;

% --- 1. �޸��ļ�·���������ڸ���3�͸���4 (��) ---
file1.path = '����3.xlsx';
file1.angle = 10;
file2.path = '����4.xlsx';
file2.angle = 15;

% --- 2. �������ṩ��ͼƬ��Ϣ���������ò��� ---
config.waveN_fit_min = 1000;      % ��������ڽϵͲ�����Ȼ���������ʵ��ſ�
% [cite: 1]
config.n1_init = 3.43;            % ����FFT����Ĺ����Ӳ�ƽ�������� (ȡ1.5-5��m��Χ���Ƽ�ֵ)
% [cite: 1]
config.n2_real_init = 3.43;         % �ĵ�������ʵ�� n2 ��ϳ�ֵ
config.k2_imag_init = 10;         % ��ĵ�������Խ�������һ����С�ĳ�ֵ
% [cite: 1]
% ����ͼƬ�ṩ�Ĺ�Sellmeier��ʽ n^2-1 = B1*L^2/(L^2-C1_sqrt^2) + ... ���ó�ֵ
% ���� B��Ӧ����ϵ��, C��Ӧ��ĸ������ƽ�� (C = lambda^2)
B1 = 10.6684; B2 = 0.00304; B3 = 1.5413; % [cite: 1]
C1_sqrt = 0.3015; C2_sqrt = 1.1348; C3_sqrt = 1104; % [cite: 1]
config.selParam_init = [B1, B2, B3, C1_sqrt^2, C2_sqrt^2, C3_sqrt^2];


disp('�����ȳ�ֵ');
data1 = readmatrix(file1.path);
thk_init1 = fft_thk_estimate(data1(:,1), data1(:,2)/100, config.n1_init, config.waveN_fit_min, file1.angle);

data2 = readmatrix(file2.path);
thk_init2 = fft_thk_estimate(data2(:,1), data2(:,2)/100, config.n1_init, config.waveN_fit_min, file2.angle);

config.thk_init = mean([thk_init1, thk_init2]);
fprintf('�ļ�3 (10��) FFT������: %.2f ��m\n', thk_init1);
fprintf('�ļ�4 (15��) FFT������: %.2f ��m\n', thk_init2);
fprintf('ƽ����ȳ�ֵ: %.2f ��m\n', config.thk_init);

fprintf('\n');
disp('--- ��ʼ������3 (10��) ---')
[result1.thk, result1.n2_complex, result1.params] = process_file(data1, file1.angle, config);

fprintf('\n');
disp('--- ��ʼ������4 (15��) ---')
[result2.thk, result2.n2_complex, result2.params] = process_file(data2, file2.angle, config);

analyze_results(result1, result2, file1.angle, file2.angle);

%% ���Ĵ�����
function [thk, n2_complex, selParam] = process_file(data, incident_angle, config)
    waveNum = data(:, 1);  % ȫ���β���
    R = data(:, 2) / 100;  % ȫ���η�����
    waveLen_full = 10000 ./ waveNum;

    filter = waveNum > config.waveN_fit_min;
    waveNum_fit = waveNum(filter);
    R_fit = R(filter);
    waveLen_fit = waveLen_full(filter);

    % ����˳��: [���, n2ʵ��, k2�鲿, 6��Sellmeierϵ��]
    x0 = [config.thk_init, config.n2_real_init, config.k2_imag_init, config.selParam_init];
    
    % �����߽� (Ϊ����ϵ���)
    lb = [config.thk_init*0.8, 3.0, 10,    10, 0.001, 1,   0.09, 1.2, 1e5]; % Lower bounds
    ub = [config.thk_init*1.2, 4.0, 30,  11, 0.005, 2,   0.1,  1.3, 2e6]; % Upper bounds

    [x_optimal, R_squared_fit] = global_fit_all_parameters(x0, lb, ub, waveNum_fit, R_fit, waveLen_fit, incident_angle);

    thk = x_optimal(1);
    n2_complex = x_optimal(2) + 1i * x_optimal(3); % ��ϳɸ���������
    selParam = x_optimal(4:9);

    fprintf('  ���պ��: %.2f ��m\n', thk);
    fprintf('  ���ճĵ׸������� ?2: %.3f + %.4fi\n', real(n2_complex), imag(n2_complex));
    fprintf('  ��������������Ŷ� R?: %.4f\n', R_squared_fit);

    n1_all = calculate_sellmeier_n(waveLen_full, selParam);
    plot_fitting_results(waveNum, R, n1_all, thk, n2_complex, incident_angle, config.waveN_fit_min);
    plot_refractive_index(waveLen_full, waveNum, n1_all);
end

%% ȫ���Ż�����
function [x_optimal, R_squared] = global_fit_all_parameters(x0, lb, ub, waveNum_fit, R_fit, waveLen_fit, incident_angle)
    theta0_rad = incident_angle * pi / 180;
    model_func = @(x, k) model_R(x, k, waveLen_fit, theta0_rad);

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
    if nargin < 5, theta0_deg = 0; end
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
    r_windowed = r_uniform(:) .* window(:);
    fft_result = fft(r_windowed);
    fft_power = abs(fft_result(1:N/2)).^2;
    dk = mean(diff(k_uniform));
    thk_axis = (0:N/2-1) * 10000 / (2 * n_avg * cos_theta1 * N * dk);
    search_range = find(thk_axis > 5 & thk_axis < 200);
    if isempty(search_range), thk = 20; return; end
    [~, max_idx_in_range] = max(fft_power(search_range));
    max_idx_rough = search_range(max_idx_in_range);
    correctNum = 3; DatePower1 = 0; DatePower2 = 0;
    for i = -correctNum:correctNum
        idx = max_idx_rough + i;
        if idx >= 1 && idx <= length(fft_power)
            power = fft_power(idx);
            DatePower1 = DatePower1 + idx * power;
            DatePower2 = DatePower2 + power;
        end
    end
    if DatePower2 > 0, f_corrected = DatePower1 / DatePower2;
    else, f_corrected = max_idx_rough; end
    thk = (f_corrected - 1) * 10000 / (2 * n_avg * cos_theta1 * N * dk);
end

%% �����������
function analyze_results(res1, res2, angle1, angle2)
    disp('--- ���ս���Ա� ---');
    fprintf('�ļ�3 (%d��) -> ���: %.2f ��m, �ĵ� ?2: %.3f + %.4fi\n', angle1, res1.thk, real(res1.n2_complex), imag(res1.n2_complex));
    fprintf('�ļ�4 (%d��) -> ���: %.2f ��m, �ĵ� ?2: %.3f + %.4fi\n', angle2, res2.thk, real(res2.n2_complex), imag(res2.n2_complex));

    thk1 = res1.thk; thk2 = res2.thk;
    thk_diff_percent = abs(thk1-thk2)/mean([thk1,thk2])*100;
    fprintf('�����ǶȲ�õĺ�ȷֱ�Ϊ %.2f ��m �� %.2f ��m����Բ���Ϊ %.2f%%��\n', thk1, thk2, thk_diff_percent);
    
    if thk_diff_percent < 2
        fprintf('  -> ����: ���һ�������ã�����ɿ���\n');
        fprintf('  -> �Ƽ����ֵ: %.2f �� %.2f ��m\n', mean([thk1, thk2]), std([thk1, thk2]));
    else
        fprintf('  -> ����: ��Ȳ���ϴ󣬽�����ģ�ͻ����ݡ�\n');
    end
end

%% --- ����ģ�����ͼ���� ---

function R = model_R(params, ~, waveLen, theta0) 
    thk = params(1);
    n2_real = params(2);
    k2_imag = params(3);
    selParam = params(4:9);
    n1 = calculate_sellmeier_n(waveLen, selParam);
    n2_complex = n2_real + 1i * k2_imag; 
    waveNum_model = 10000 ./ waveLen;
    R = model_R_vectorized(thk, n1, n2_complex, waveNum_model, theta0);
end

function R = model_R_vectorized(thk, n1, n2_complex, waveNum, theta0)
    n0 = 1.0; 
    n1 = n1(:);
    n2 = n2_complex(:);
    waveNum = waveNum(:);

    sin_theta1 = n0 * sin(theta0) ./ n1;
    cos_theta1 = sqrt(1 - sin_theta1.^2);
    sin_theta2 = n0 * sin(theta0) ./ n2;
    cos_theta2 = sqrt(1 - sin_theta2.^2);

    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    r12_s = (n1.*cos_theta1 - n2.*cos_theta2) ./ (n1.*cos_theta1 + n2.*cos_theta2);
    r01_p = (n1*cos(theta0) - n0.*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    r12_p = (n2.*cos_theta1 - n1.*cos_theta2) ./ (n2.*cos_theta1 + n1.*cos_theta2);

    delta = 4 * pi * n1 .* thk .* cos_theta1 .* waveNum / 10000;
    
    r_s_total = (r01_s + r12_s .* exp(1i*delta)) ./ (1 + r01_s .* r12_s .* exp(1i*delta));
    r_p_total = (r01_p + r12_p .* exp(1i*delta)) ./ (1 + r01_p .* r12_p .* exp(1i*delta));

    R_s = abs(r_s_total).^2;
    R_p = abs(r_p_total).^2;

    R = (R_s + R_p) / 2;
    R = real(R(:));
end

function n = calculate_sellmeier_n(waveLen, params)
    B1 = params(1); B2 = params(2); B3 = params(3);
    C1 = params(4); C2 = params(5); C3 = params(6);
    lambda_sq = waveLen.^2; 
    eps = 1e-10;
    term1 = B1 * lambda_sq ./ (lambda_sq - C1 + eps);
    term2 = B2 * lambda_sq ./ (lambda_sq - C2 + eps);
    term3 = B3 * lambda_sq ./ (lambda_sq - C3 + eps);
    n_squared = 1 + term1 + term2 + term3;
    n_squared(n_squared < 1) = 1;
    n = real(sqrt(n_squared)); 
    n(n < 1) = 1;
end

function plot_refractive_index(waveLen, waveNum, n1_all)
    figure('Name', '��ϵõ������Ӳ�������ɫɢ����', 'Position', [100, 100, 1200, 500]);
    subplot(1,2,1); plot(waveLen, n1_all, 'b-', 'LineWidth', 2); xlabel('���� (��m)'); ylabel('���Ӳ������� n_1'); title('������ vs ����'); grid on;
    subplot(1,2,2); plot(waveNum, n1_all, 'r-', 'LineWidth', 2); xlabel('���� (cm^{-1})'); ylabel('���Ӳ������� n_1'); title('������ vs ����'); grid on; xlim([min(waveNum), max(waveNum)]);
end

function plot_fitting_results(waveNum, R, n1, thk, n2_complex, incident_angle, waveN_fit_min)
    theta0 = incident_angle * pi / 180;
    R_fitted = model_R_vectorized(thk, n1, n2_complex, waveNum, theta0);
    R_squared_full = 1 - sum((R - R_fitted).^2) / sum((R - mean(R)).^2);

    figure('Name', ['��Ͻ������ (����� ' num2str(incident_angle) '��)'], 'Position', [150, 150, 1200, 600]);
    plot(waveNum, R*100, 'b.', 'MarkerSize', 5, 'DisplayName', 'ʵ������');
    hold on;
    plot(waveNum, R_fitted*100, 'r-', 'LineWidth', 2, 'DisplayName', '���ģ��');
    ylim_vals = ylim;
    h = fill([waveN_fit_min, max(waveNum), max(waveNum), waveN_fit_min], [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], 'k', 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'DisplayName', '�������');
    uistack(h, 'bottom');
    xlabel('���� (cm^{-1})'); ylabel('������ (%)');
    title_str = sprintf('��������Ͻ�� (��=%d��, d=%.2f ��m, ?_2=%.3f+%.4fi, ȫ��R?=%.4f)', ...
        incident_angle, thk, real(n2_complex), imag(n2_complex), R_squared_full);
    title(title_str);
    legend('Location', 'best'); grid on; xlim([min(waveNum), max(waveNum)]);
end