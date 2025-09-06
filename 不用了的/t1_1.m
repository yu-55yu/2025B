%% ��������̼�������Ӳ��Ȳ���
% �ļ�����main_sic_thickness.m
function main_sic_thickness()
    clc; clear; close all;
    
    %% ������1��10������ǣ�
    fprintf('==================================================\n');
    fprintf('������1���ݣ������10�㣩\n');
    fprintf('==================================================\n');
    
    [thickness1, n2_1, params1, rmse1] = process_data('����1.xlsx', 10);
    
    fprintf('\n����1���:\n');
    fprintf('���Ӳ���: %.2f ��m\n', thickness1);
    fprintf('�ĵ�������: %.3f\n', n2_1);
    fprintf('������Ү������:\n');
    fprintf('  B1=%.3f, B2=%.3f, B3=%.3f\n', params1(1), params1(2), params1(3));
    fprintf('  C1=%.3f, C2=%.3f, C3=%.3f\n', params1(4), params1(5), params1(6));
    fprintf('RMSE: %.4f\n', rmse1);
    
    %% ������2��15������ǣ�
    fprintf('\n==================================================\n');
    fprintf('������2���ݣ������15�㣩\n');
    fprintf('==================================================\n');
    
    [thickness2, n2_2, params2, rmse2] = process_data('����2.xlsx', 15);
    
    fprintf('\n����2���:\n');
    fprintf('���Ӳ���: %.2f ��m\n', thickness2);
    fprintf('�ĵ�������: %.3f\n', n2_2);
    fprintf('������Ү������:\n');
    fprintf('  B1=%.3f, B2=%.3f, B3=%.3f\n', params2(1), params2(2), params2(3));
    fprintf('  C1=%.3f, C2=%.3f, C3=%.3f\n', params2(4), params2(5), params2(6));
    fprintf('RMSE: %.4f\n', rmse2);
    
    %% �ɿ��Է���
    fprintf('\n==================================================\n');
    fprintf('�ɿ��Է���:\n');
    fprintf('==================================================\n');
    
    thickness_diff = abs(thickness1 - thickness2);
    relative_diff = thickness_diff / ((thickness1 + thickness2) / 2) * 100;
    
    fprintf('���ֽǶȲ�õĺ�Ȳ�: %.2f ��m\n', thickness_diff);
    fprintf('���ƫ��: %.1f%%\n', relative_diff);
    
    if thickness_diff < 1.0
        fprintf('����ɿ���: ���㣨��Ȳ�<1��m��\n');
    elseif thickness_diff < 2.0
        fprintf('����ɿ���: ���ã���Ȳ�<2��m��\n');
    else
        fprintf('����ɿ���: ��Ҫ��һ����֤\n');
    end
end

%% �����������ļ�
function [thickness, n2, sellmeier_params, rmse] = process_data(filename, incident_angle)
    % ��ȡ����
    data = readmatrix(filename);
    wavenumber = data(:, 1);  % ���� cm^-1
    reflectance = data(:, 2) / 100;  % ������ת��ΪС��
    
    % ����(��m) = 10000/����
    wavelength = 10000 ./ wavenumber;
    
    % �����ת��Ϊ����
    theta0 = incident_angle * pi / 180;
    
    % ʹ��FFT���ƺ�ȳ�ֵ
    thickness_init = estimate_thickness_fft(wavenumber, reflectance);
    fprintf('FFT���ƺ�ȳ�ֵ: %.2f ��m\n', thickness_init);
    
    % �����Ż�����
    % ����˳��: [B1, B2, B3, C1, C2, C3, thickness, n2]
    x0 = [2.0, 0.5, 0.1, 0.01, 0.1, 10, thickness_init, 2.55];
    lb = [0.1, 0.01, 0.001, 0.001, 0.01, 1, 1, 2.0];
    ub = [10, 5, 2, 1, 10, 100, 200, 3.0];
    
    % ����Ŀ�꺯��
    objective = @(x) sum(calculate_residuals(x, wavelength, wavenumber, reflectance, theta0).^2);
    
    % ȫ���Ż�ѡ��
    fprintf('��ʼȫ���Ż�...\n');
    
    % ʹ���Ŵ��㷨����ȫ���Ż�
    options_ga = optimoptions('ga', ...
        'PopulationSize', 100, ...
        'MaxGenerations', 50, ...
        'Display', 'iter', ...
        'UseParallel', true, ...
        'FunctionTolerance', 1e-8);
    
    [x_ga, ~] = ga(objective, 8, [], [], [], [], lb, ub, [], options_ga);
    
    % ʹ�þֲ��Ż���һ����ϸ��
    fprintf('\n��ʼ�ֲ��Ż�...\n');
    options_local = optimoptions('fmincon', ...
        'Display', 'iter', ...
        'Algorithm', 'interior-point', ...
        'MaxIterations', 500, ...
        'OptimalityTolerance', 1e-10, ...
        'StepTolerance', 1e-10);
    
    [x_optimal, fval] = fmincon(objective, x_ga, [], [], [], [], lb, ub, [], options_local);
    
    % ��ȡ���
    sellmeier_params = x_optimal(1:6);
    thickness = x_optimal(7);
    n2 = x_optimal(8);
    
    % ����RMSE
    residuals = calculate_residuals(x_optimal, wavelength, wavenumber, reflectance, theta0);
    rmse = sqrt(mean(residuals.^2));
    
    % ���ƽ��
    plot_results(wavenumber, reflectance, x_optimal, wavelength, theta0, incident_angle);
end

%% FFT���ƺ��
function thickness_est = estimate_thickness_fft(wavenumber, reflectance)
    % ���Ļ�����
    reflectance_centered = reflectance - mean(reflectance);
    
    % FFT����
    N = length(wavenumber);
    fft_result = fft(reflectance_centered);
    
    % ����Ƶ��
    delta_k = mean(diff(wavenumber));
    freq = (0:N-1) / (N * delta_k);
    
    % �ҵ���Ƶ�ʣ��ų�DC������
    fft_power = abs(fft_result).^2;
    [~, max_idx] = max(fft_power(2:floor(N/2)));
    main_freq = freq(max_idx + 1);
    
    % ���ƺ�� d �� 10000/(2n*��Ƶ��)
    n_avg = 2.6;  % ̼�������������
    thickness_est = 10000 / (2 * n_avg * main_freq);
    
    % �����ں���Χ��
    thickness_est = max(10, min(150, thickness_est));
end

%% ����������Ү��������
function n = sellmeier(wavelength, params)
    B1 = params(1); B2 = params(2); B3 = params(3);
    C1 = params(4); C2 = params(5); C3 = params(6);
    
    lambda_sq = wavelength.^2;
    
    % �������
    term1 = B1 * lambda_sq ./ (lambda_sq - C1 + 1e-10);
    term2 = B2 * lambda_sq ./ (lambda_sq - C2 + 1e-10);
    term3 = B3 * lambda_sq ./ (lambda_sq - C3 + 1e-10);
    
    n_squared = 1 + term1 + term2 + term3;
    n = sqrt(abs(n_squared));
end

%% ���㷴����ģ��
function R = reflectance_model(params, wavelength, wavenumber, theta0)
    % ��ȡ����
    sellmeier_params = params(1:6);
    thickness = params(7);
    n2 = params(8);
    
    % �������Ӳ�������
    n1 = sellmeier(wavelength, sellmeier_params);
    n0 = 1.0;  % ����������
    
    % ��������ǣ�˹�������ɣ�
    sin_theta1 = n0 * sin(theta0) ./ n1;
    sin_theta2 = n0 * sin(theta0) / n2;
    
    % ȷ��ֵ����Ч��Χ��
    sin_theta1 = max(-1, min(1, sin_theta1));
    sin_theta2 = max(-1, min(1, sin_theta2));
    
    cos_theta0 = cos(theta0);
    cos_theta1 = sqrt(1 - sin_theta1.^2);
    cos_theta2 = sqrt(1 - sin_theta2^2);
    
    % ���������ϵ����sƫ��
    r01 = (n0*cos_theta0 - n1.*cos_theta1) ./ (n0*cos_theta0 + n1.*cos_theta1);
    t01 = 2*n0*cos_theta0 ./ (n0*cos_theta0 + n1.*cos_theta1);
    t10 = 2*n1.*cos_theta1 ./ (n0*cos_theta0 + n1.*cos_theta1);
    r12 = (n1.*cos_theta1 - n2*cos_theta2) ./ (n1.*cos_theta1 + n2*cos_theta2);
    
    % ��λ�� �� = 4��n�6�9d��cos(�ȁ6�9)����/10000
    delta = 4 * pi * n1 .* thickness .* cos_theta1 .* wavenumber / 10000;
    
    % ���㷴����
    numerator = r01.^2 + r12^2 * t01.^2 .* t10.^2 + 2*r01.*r12.*t01.*t10.*cos(delta);
    denominator = 1 + r01.^2 * r12^2 + 2*r01.*r12.*cos(delta);
    
    R = numerator ./ denominator;
end

%% ����в�
function residuals = calculate_residuals(params, wavelength, wavenumber, reflectance, theta0)
    R_model = reflectance_model(params, wavelength, wavenumber, theta0);
    % ��Ȩ�в�
    weights = sqrt(reflectance + 0.01);
    residuals = (R_model - reflectance) .* weights;
end

%% ���ƽ��
function plot_results(wavenumber, reflectance, params, wavelength, theta0, incident_angle)
    % ������ϵķ�����
    R_fitted = reflectance_model(params, wavelength, wavenumber, theta0);
    
    % ��ȡ����
    thickness = params(7);
    n2 = params(8);
    
    % ����ͼ��
    figure('Position', [100, 100, 1000, 600]);
    
    % ��ͼ1�������ʶԱ�
    subplot(2, 1, 1);
    plot(wavenumber, reflectance, 'b.', 'MarkerSize', 2);
    hold on;
    plot(wavenumber, R_fitted, 'r-', 'LineWidth', 1.5);
    xlabel('���� (cm^{-1})');
    ylabel('������');
    title(sprintf('��������Ͻ�� (�����=%d��, ���=%.2f ��m, n_2=%.3f)', ...
        incident_angle, thickness, n2));
    legend('ʵ������', '��Ͻ��', 'Location', 'best');
    grid on;
    grid minor;
    
    % ��ͼ2���в�ֲ�
    subplot(2, 1, 2);
    residuals = reflectance - R_fitted;
    plot(wavenumber, residuals, 'g.', 'MarkerSize', 2);
    hold on;
    yline(0, 'k--', 'LineWidth', 1);
    xlabel('���� (cm^{-1})');
    ylabel('�в�');
    title(sprintf('�в�ֲ� (RMSE = %.4f)', sqrt(mean(residuals.^2))));
    grid on;
    grid minor;
    
    % ��������
    sgtitle(sprintf('̼�������Ӳ��Ȳ������ (����� %d��)', incident_angle));
end

%% ��ȡ���в�����Ӧ��������
function n1_all = get_n1_all(wavelength, sellmeier_params)
    n1_all = sellmeier(wavelength, sellmeier_params);
end

%% ����ĸ�����������������ͷ���
function batch_analysis()
    % �˺������������������ļ�����в��������Է���
    
    % �Ƕ�ɨ�����
    angles = 5:5:30;
    thicknesses = zeros(length(angles), 1);
    
    for i = 1:length(angles)
        % ��������ж�Ӧ�Ƕȵ������ļ�
        filename = sprintf('data_angle_%d.xlsx', angles(i));
        if exist(filename, 'file')
            [thickness, ~, ~, ~] = process_data(filename, angles(i));
            thicknesses(i) = thickness;
        end
    end
    
    % ���ƽǶ�������
    figure;
    plot(angles, thicknesses, 'bo-', 'LineWidth', 2);
    xlabel('����� (��)');
    ylabel('������� (��m)');
    title('��Ȳ����ĽǶ�������');
    grid on;
end