%% ������������Ү��������Ϻ����Ӳ��ȼ���
% �ļ�����sellmeier_thickness_solver_improved.m

function sellmeier_thickness_solver_improved()
    clc; clear; close all;
    
    %% ������1��10������ǣ�
    fprintf('==================================================\n');
    fprintf('������1���ݣ������10�㣩\n');
    fprintf('==================================================\n');
    
    [thickness1, n1_all_1, params1] = process_file('����1.xlsx', 10);
    
    %% ������2��15������ǣ�  
    fprintf('\n==================================================\n');
    fprintf('������2���ݣ������15�㣩\n');
    fprintf('==================================================\n');
    
    [thickness2, n1_all_2, params2] = process_file('����2.xlsx', 15);
    
    %% ����Ա�����֤
    fprintf('\n==================================================\n');
    fprintf('����Աȷ���\n');
    fprintf('==================================================\n');
    fprintf('����1���: %.2f ��m\n', thickness1);
    fprintf('����2���: %.2f ��m\n', thickness2);
    fprintf('��Ȳ���: %.2f ��m (%.1f%%)\n', ...
        abs(thickness1-thickness2), ...
        abs(thickness1-thickness2)/mean([thickness1,thickness2])*100);
    
    % ����ɿ��Է���
    analyze_reliability(thickness1, thickness2, params1, params2);
end

%% �������ļ�
function [thickness, n1_all, sellmeier_params] = process_file(filename, incident_angle)
    
    %% 1. ��ȡ����
    data = readmatrix(filename);
    wavenumber = data(:, 1);  % ���� cm^-1
    reflectance = data(:, 2) / 100;  % ������ת��ΪС��
    wavelength = 10000 ./ wavenumber;  % ���� ��m
    
    %% 2. ��һ�����ø߲����������������Ү������
    fprintf('\n����1: ʹ�ø߲�������(>2000 cm^-1)���������Ү������...\n');
    
    % ѡ��߲������ݣ�ֻ������ϲ�����
    high_k_mask = wavenumber > 2000;
    high_k_wavelength = wavelength(high_k_mask);
    high_k_reflectance = reflectance(high_k_mask);
    
    % ���������Ү���������Ľ��棩
    sellmeier_params = fit_sellmeier_params_improved(high_k_wavelength, high_k_reflectance, incident_angle);
    
    fprintf('������Ү������:\n');
    fprintf('  B1=%.4f, B2=%.4f, B3=%.4f\n', sellmeier_params(1), sellmeier_params(2), sellmeier_params(3));
    fprintf('  C1=%.4f, C2=%.4f, C3=%.4f\n', sellmeier_params(4), sellmeier_params(5), sellmeier_params(6));
    
    %% 3. �������в�����������
    n1_all = calculate_sellmeier_n(wavelength, sellmeier_params);
    
    % ��ʾ������ͳ��
    display_n_statistics(n1_all, wavenumber);
    
    % �������������ߣ�ȫ���Σ�
    plot_refractive_index(wavelength, wavenumber, n1_all, high_k_mask);
    
    %% 4. �ڶ�����ʹ��FFT���ƺ�ȳ�ֵ
    fprintf('\n����2: FFT���ƺ�ȳ�ֵ...\n');
    
    % ʹ��ƽ��������
    n_avg = mean(n1_all);
    
    % FFT�������Ľ��棩
    thickness_init = fft_thickness_estimate_improved(wavenumber, reflectance, n_avg);
    fprintf('FFT���ƺ��: %.2f ��m\n', thickness_init);
    
    %% 5. ����������ȷ��Ϻ��
    fprintf('\n����3: ��ȷ��Ϻ��...\n');
    
    % �ĵ������ʳ�ֵ��̼�������ֵ��
    n2_init = 2.55;
    
    % �Ż���Ⱥͳĵ������ʣ��Ľ��棩
    [thickness, n2, fit_quality] = fit_thickness_improved(wavenumber, reflectance, n1_all, incident_angle, thickness_init, n2_init);
    
    fprintf('���պ��: %.2f ��m\n', thickness);
    fprintf('�ĵ�������: %.3f\n', n2);
    fprintf('����Ŷ� R?: %.4f\n', fit_quality);
    
    %% 6. ������Ͻ��
    plot_fitting_results_improved(wavenumber, reflectance, n1_all, thickness, n2, incident_angle);
end

%% �Ľ���������Ү���������
function params = fit_sellmeier_params_improved(wavelength, reflectance, incident_angle)
    
    % �����ʼֵ���ԣ�ѡ�����Ž��
    initial_guesses = [
        [5.0, 0.5, 0.1, 0.01, 0.5, 20];    % �����²�
        [6.5, 0.3, 0.05, 0.02, 0.3, 15];   % ̼�������ֵ
        [4.0, 0.8, 0.2, 0.005, 0.8, 25];   % ����²�
    ];
    
    % �����߽磨�����ɣ�
    lb = [0.1, 0.001, 0.0001, 0.0001, 0.001, 0.1];
    ub = [20, 10, 5, 2, 20, 200];
    
    theta0 = incident_angle * pi / 180;
    
    best_params = [];
    best_resnorm = inf;
    
    % ���Բ�ͬ��ʼֵ
    for i = 1:size(initial_guesses, 1)
        x0 = initial_guesses(i, :);
        
        objective = @(x) sellmeier_objective_improved(x, wavelength, reflectance, theta0);
        
        options = optimoptions('lsqnonlin', ...
            'Display', 'off', ...
            'MaxIterations', 2000, ...
            'FunctionTolerance', 1e-12, ...
            'StepTolerance', 1e-12, ...
            'Algorithm', 'trust-region-reflective');
        
        try
            [params_temp, resnorm_temp] = lsqnonlin(objective, x0, lb, ub, options);
            
            if resnorm_temp < best_resnorm
                best_params = params_temp;
                best_resnorm = resnorm_temp;
            end
        catch
            continue;
        end
    end
    
    params = best_params;
    fprintf('������ϲв�: %.6f\n', best_resnorm);
end

%% �Ľ���������Ү��Ŀ�꺯��
function residuals = sellmeier_objective_improved(params, wavelength, reflectance, theta0)
    
    % ����������
    n1 = calculate_sellmeier_n(wavelength, params);
    
    % ȷ���������ں���Χ
    if any(n1 < 1.5) || any(n1 > 4.0) || any(imag(n1) ~= 0)
        residuals = ones(size(reflectance)) * 1e6;  % �ͷ���
        return;
    end
    
    % �����ͳĵ׵�������
    n0 = 1.0;
    
    % ��������ǣ�˹�������ɣ�
    sin_theta1 = n0 * sin(theta0) ./ n1;
    cos_theta1 = real(sqrt(1 - sin_theta1.^2));
    
    % ��ƫ���ķ���������ϵ��
    % sƫ��
    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    
    % pƫ��
    r01_p = (n1*cos(theta0) - n0.*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    
    % ��ƫ��ⷴ����
    R_surface = (abs(r01_s).^2 + abs(r01_p).^2) / 2;
    
    % �вʹ���������Ȩ��
    weight = 1 ./ (reflectance + 0.01);  % �������
    residuals = (R_surface - reflectance) .* sqrt(weight);
end

%% �Ľ���FFT��ȹ���
function thickness = fft_thickness_estimate_improved(wavenumber, reflectance, n_avg)
    
    % ����Ԥ����
    reflectance_smooth = smooth(reflectance, 5);  % ƽ������
    reflectance_ac = reflectance_smooth - mean(reflectance_smooth);
    
    % ��ֵ�����ܵľ�������
    N = 2^nextpow2(4*length(wavenumber));  % ʹ��2���ݴ������FFTЧ��
    k_uniform = linspace(min(wavenumber), max(wavenumber), N);
    r_uniform = interp1(wavenumber, reflectance_ac, k_uniform, 'pchip');
    
    % �Ӵ�����Ƶ��й¶
    window = hann(N);
    r_windowed = r_uniform(:) .* window(:);
    
    % FFT����
    fft_result = fft(r_windowed);
    fft_power = abs(fft_result).^2;
    
    % Ƶ���ᣨ�����
    dk = mean(diff(k_uniform));
    thickness_axis = (0:N-1) * 10000 / (2 * n_avg * N * dk);
    
    % �ҵ����壨�ų���Ƶ��������
    search_range = find(thickness_axis > 5 & thickness_axis < 200);
    [~, max_idx] = max(fft_power(search_range));
    thickness = thickness_axis(search_range(max_idx));
end

%% �Ľ��ĺ�����
function [thickness, n2, R_squared] = fit_thickness_improved(wavenumber, reflectance, n1, incident_angle, d_init, n2_init)
    
    theta0 = incident_angle * pi / 180;
    
    % ����Ŀ�꺯��
    model_func = @(x, xdata) model_reflectance_vectorized(x, xdata, n1, theta0);
    
    % ��ʼֵ�ͱ߽�
    x0 = [d_init, n2_init];
    lb = [0.5, 2.0];   % ����ķ�Χ
    ub = [300, 3.5];
    
    % �Ż�ѡ��
    options = optimoptions('lsqcurvefit', ...
        'Display', 'iter', ...
        'MaxIterations', 1000, ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'Algorithm', 'trust-region-reflective', ...
        'UseParallel', true);  % ���м���
    
    % ���
    [x_optimal, resnorm] = lsqcurvefit(model_func, x0, wavenumber, reflectance, lb, ub, options);
    
    thickness = x_optimal(1);
    n2 = x_optimal(2);
    
    % ��������Ŷ�
    R_fitted = model_func(x_optimal, wavenumber);
    SS_tot = sum((reflectance - mean(reflectance)).^2);
    SS_res = sum((reflectance - R_fitted).^2);
    R_squared = 1 - SS_res/SS_tot;
    
    fprintf('������ϲв�: %.6f\n', resnorm);
end

%% �������ķ�����ģ��
function R = model_reflectance_vectorized(params, wavenumber, n1, theta0)
    
    thickness = params(1);
    n2 = params(2);
    n0 = 1.0;
    
    % ȷ��������������
    n1 = n1(:);
    wavenumber = wavenumber(:);
    
    % ������������
    sin_theta1 = n0 * sin(theta0) ./ n1;
    cos_theta1 = real(sqrt(1 - sin_theta1.^2));
    
    sin_theta2 = n0 * sin(theta0) / n2;
    cos_theta2 = real(sqrt(1 - sin_theta2^2));
    
    % ��ƫ�����ۺϷ���ϵ��
    % sƫ��
    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    r12_s = (n1.*cos_theta1 - n2*cos_theta2) ./ (n1.*cos_theta1 + n2*cos_theta2);
    
    % pƫ��
    r01_p = (n1*cos(theta0) - n0.*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    r12_p = (n2*cos_theta1 - n1.*cos_theta2) ./ (n2*cos_theta1 + n1.*cos_theta2);
    
    % ��λ��
    delta = 4 * pi * n1 .* thickness .* cos_theta1 .* wavenumber / 10000;
    
    % ˫�������湫ʽ��s��p��ƽ����
    R_s = abs((r01_s + r12_s .* exp(1i*delta)) ./ (1 + r01_s .* r12_s .* exp(1i*delta))).^2;
    R_p = abs((r01_p + r12_p .* exp(1i*delta)) ./ (1 + r01_p .* r12_p .* exp(1i*delta))).^2;
    
    R = (R_s + R_p) / 2;
    R = real(R(:));  % ȷ�������ʵ��������
end

%% ������Ү�����̼���������
function n = calculate_sellmeier_n(wavelength, params)
    % Sellmeier����: n^2 = 1 + ��[Bi*��^2/(��^2-Ci)]
    B1 = params(1); B2 = params(2); B3 = params(3);
    C1 = params(4); C2 = params(5); C3 = params(6);
    
    lambda_sq = wavelength.^2;
    
    % �������Ľ���ֵ�ȶ��ԣ�
    eps = 1e-10;  % �������
    term1 = B1 * lambda_sq ./ (lambda_sq - C1 + eps);
    term2 = B2 * lambda_sq ./ (lambda_sq - C2 + eps);  
    term3 = B3 * lambda_sq ./ (lambda_sq - C3 + eps);
    
    % ȷ��ÿ�����ֵ����
    term1 = max(0, real(term1));
    term2 = max(0, real(term2));
    term3 = max(0, real(term3));
    
    n_squared = 1 + term1 + term2 + term3;
    
    % ȷ��nΪʵ���Һ���
    n = real(sqrt(n_squared));
    n = max(1.0, min(4.0, n));  % �����ں���Χ
end

%% ��������������
function plot_refractive_index(wavelength, wavenumber, n1_all, high_k_mask)
    
    figure('Position', [100, 100, 1200, 600]);
    
    % ��ͼ1��������vs������ȫ���Σ�
    subplot(1,2,1);
    plot(wavelength, n1_all, 'b-', 'LineWidth', 2);
    hold on;
    plot(wavelength(high_k_mask), n1_all(high_k_mask), 'r.', 'MarkerSize', 8);
    xlabel('���� (��m)');
    ylabel('������ n');
    title('���Ӳ�������ɫɢ����');
    legend('ȫ����������', '�߲����������', 'Location', 'best');
    grid on;
    ylim([min(n1_all)*0.98, max(n1_all)*1.02]);
    
    % ��ͼ2��������vs������ȫ���Σ�
    subplot(1,2,2);
    plot(wavenumber, n1_all, 'b-', 'LineWidth', 2);
    hold on;
    plot(wavenumber(high_k_mask), n1_all(high_k_mask), 'r.', 'MarkerSize', 8);
    xline(2000, 'k--', '��ϱ߽�', 'LineWidth', 1.5);
    xlabel('���� (cm^{-1})');
    ylabel('������ n');
    title('���Ӳ������ʣ�������');
    legend('ȫ����������', '�߲����������', 'Location', 'best');
    grid on;
    xlim([min(wavenumber), max(wavenumber)]);
    
    sgtitle('������Ү��������Ͻ��');
end

%% �Ľ�����Ͻ����ͼ
function plot_fitting_results_improved(wavenumber, reflectance, n1, thickness, n2, incident_angle)
    
    theta0 = incident_angle * pi / 180;
    
    % ������ϵķ�����
    R_fitted = model_reflectance_vectorized([thickness, n2], wavenumber, n1, theta0);
    
    % �����ۺ�ͼ��
    figure('Position', [100, 100, 1400, 900]);
    
    % ��ͼ1��ȫ���η����ʶԱ�
    subplot(3,2,[1,2]);
    plot(wavenumber, reflectance*100, 'b.', 'MarkerSize', 2);
    hold on;
    plot(wavenumber, R_fitted*100, 'r-', 'LineWidth', 1.5);
    xlabel('���� (cm^{-1})');
    ylabel('������ (%)');
    title(sprintf('��������Ͻ�� (��=%d��, d=%.2f ��m, n?=%.3f)', ...
        incident_angle, thickness, n2));
    legend('ʵ������', '��Ͻ��', 'Location', 'best');
    grid on;
    xlim([min(wavenumber), max(wavenumber)]);
    
    % ��ͼ2���в����
    subplot(3,2,3);
    residuals = (reflectance - R_fitted) * 100;
    plot(wavenumber, residuals, 'g.', 'MarkerSize', 2);
    hold on;
    yline(0, 'k--', 'LineWidth', 1);
    yline(std(residuals), 'r--', '+��');
    yline(-std(residuals), 'r--', '-��');
    xlabel('���� (cm^{-1})');
    ylabel('�в� (%)');
    title(sprintf('�в�ֲ� (RMSE = %.4f%%)', sqrt(mean(residuals.^2))));
    grid on;
    xlim([min(wavenumber), max(wavenumber)]);
    
    % ��ͼ3���в�ֱ��ͼ
    subplot(3,2,4);
    histogram(residuals, 30, 'Normalization', 'probability');
    xlabel('�в� (%)');
    ylabel('����');
    title('�в�ֲ�ֱ��ͼ');
    grid on;
    
    % ��ͼ4���߲�������ϸ��
    subplot(3,2,5);
    high_k_mask = wavenumber > 2000;
    plot(wavenumber(high_k_mask), reflectance(high_k_mask)*100, 'b.', 'MarkerSize', 3);
    hold on;
    plot(wavenumber(high_k_mask), R_fitted(high_k_mask)*100, 'r-', 'LineWidth', 1.5);
    xlabel('���� (cm^{-1})');
    ylabel('������ (%)');
    title('�߲����������ϸ�� (>2000 cm^{-1})');
    legend('ʵ������', '��Ͻ��', 'Location', 'best');
    grid on;
    
    % ��ͼ5���Ͳ�������ϸ��
    subplot(3,2,6);
    low_k_mask = wavenumber < 1000;
    plot(wavenumber(low_k_mask), reflectance(low_k_mask)*100, 'b.', 'MarkerSize', 3);
    hold on;
    plot(wavenumber(low_k_mask), R_fitted(low_k_mask)*100, 'r-', 'LineWidth', 1.5);
    xlabel('���� (cm^{-1})');
    ylabel('������ (%)');
    title('�Ͳ����������ϸ�� (<1000 cm^{-1})');
    legend('ʵ������', '��Ͻ��', 'Location', 'best');
    grid on;
    
    % ���㲢��ʾ����Ŷ�ָ��
    R_squared = 1 - sum((reflectance - R_fitted).^2) / sum((reflectance - mean(reflectance)).^2);
    MAE = mean(abs(residuals));
    
    % ����ܱ���
    sgtitle(sprintf('̼�������Ӳ��Ȳ������ (R? = %.4f, MAE = %.4f%%)', R_squared, MAE));
end

%% ��ʾ������ͳ��
function display_n_statistics(n1_all, wavenumber)
    fprintf('\n������ͳ����Ϣ:\n');
    fprintf('  ��Χ: [%.4f, %.4f]\n', min(n1_all), max(n1_all));
    fprintf('  ƽ��ֵ: %.4f\n', mean(n1_all));
    fprintf('  ��׼��: %.4f\n', std(n1_all));
    fprintf('  ��λ��: %.4f\n', median(n1_all));
    
    % �ض��������������
    specific_k = [500, 1000, 1500, 2000, 2500, 3000, 3500];
    fprintf('\n�ض��������������:\n');
    for k = specific_k
        [~, idx] = min(abs(wavenumber - k));
        if idx <= length(n1_all)
            fprintf('  %4d cm^-1: n = %.4f (�� = %.2f ��m)\n', ...
                k, n1_all(idx), 10000/k);
        end
    end
end

%% ����ɿ��Է���
function analyze_reliability(thickness1, thickness2, params1, params2)
    
    fprintf('\n=== ����ɿ��Է��� ===\n');
    
    % 1. ���һ���Լ���
    thickness_diff_percent = abs(thickness1-thickness2)/mean([thickness1,thickness2])*100;
    if thickness_diff_percent < 5
        fprintf('? ���һ�������� (���� < 5%%)\n');
    else
        fprintf('? ��Ȳ���ϴ� (%.1f%%)����Ҫ��һ�����\n', thickness_diff_percent);
    end
    
    % 2. ������Ү������������
    fprintf('\n������Ү�������Ա�:\n');
    param_names = {'B1', 'B2', 'B3', 'C1', 'C2', 'C3'};
    for i = 1:6
        diff = abs(params1(i) - params2(i))/mean([params1(i), params2(i)])*100;
        fprintf('  %s: %.4f vs %.4f (����: %.1f%%)\n', ...
            param_names{i}, params1(i), params2(i), diff);
    end
    
    % 3. ��������Լ���
    fprintf('\n��������Լ���:\n');
    
    % ����ȷ�Χ
    if thickness1 > 1 && thickness1 < 200 && thickness2 > 1 && thickness2 < 200
        fprintf('? ����ں���Χ�� (1-200 ��m)\n');
    else
        fprintf('? ��ȿ��ܳ�������Χ\n');
    end
    
    % ��������ʷ�Χ��̼�������ֵ2.5-2.8��
    lambda_test = 2.0;  % ���Բ���2��m
    n_test1 = calculate_sellmeier_n(lambda_test, params1);
    n_test2 = calculate_sellmeier_n(lambda_test, params2);
    
    if n_test1 > 2.3 && n_test1 < 3.0 && n_test2 > 2.3 && n_test2 < 3.0
        fprintf('? �������ں���Χ�� (2.3-3.0)\n');
        fprintf('  ��=2��mʱ: n1=%.3f, n2=%.3f\n', n_test1, n_test2);
    else
        fprintf('? �����ʿ����쳣\n');
    end
    
    % 4. ����
    fprintf('\n����:\n');
    if thickness_diff_percent < 5
        fprintf('? ����ɿ��������ǶȲ������һ��\n');
        fprintf('? �Ƽ����ֵ: %.2f �� %.2f ��m\n', ...
            mean([thickness1, thickness2]), ...
            abs(thickness1-thickness2)/2);
    else
        fprintf('? ������ʵ�����������²���\n');
        fprintf('? ������Ҫ���Ƕ��������ЧӦ\n');
    end
end