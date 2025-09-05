%% ������������Ү��������Ϻ����Ӳ��ȼ���
% �ļ�����sellmeier_thickness_solver.m

function sellmeier_thickness_solver()
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
    
    %% ����Ա�
    fprintf('\n==================================================\n');
    fprintf('����Աȷ���\n');
    fprintf('==================================================\n');
    fprintf('����1���: %.2f ��m\n', thickness1);
    fprintf('����2���: %.2f ��m\n', thickness2);
    fprintf('��Ȳ���: %.2f ��m (%.1f%%)\n', ...
        abs(thickness1-thickness2), ...
        abs(thickness1-thickness2)/mean([thickness1,thickness2])*100);
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
    
    % ѡ��߲�������
    high_k_mask = wavenumber > 2000;
    high_k_wavelength = wavelength(high_k_mask);
    high_k_reflectance = reflectance(high_k_mask);
    high_k_wavenumber = wavenumber(high_k_mask);
    
    % ���������Ү������
    sellmeier_params = fit_sellmeier_params(high_k_wavelength, high_k_reflectance, incident_angle);
    
    fprintf('������Ү������:\n');
    fprintf('  B1=%.4f, B2=%.4f, B3=%.4f\n', sellmeier_params(1), sellmeier_params(2), sellmeier_params(3));
    fprintf('  C1=%.4f, C2=%.4f, C3=%.4f\n', sellmeier_params(4), sellmeier_params(5), sellmeier_params(6));
    
    %% 3. �������в�����������
    n1_all = calculate_sellmeier_n(wavelength, sellmeier_params);
    
    % ��ʾ������ͳ��
    display_n_statistics(n1_all, wavenumber);
    
    % ��������������
    figure;
    subplot(2,1,1);
    plot(wavelength, n1_all, 'b-', 'LineWidth', 1.5);
    xlabel('���� (��m)');
    ylabel('������ n');
    title('���Ӳ������ʵĲ���������');
    grid on;
    
    subplot(2,1,2);
    plot(wavenumber, n1_all, 'r-', 'LineWidth', 1.5);
    xlabel('���� (cm^{-1})');
    ylabel('������ n');
    title('���Ӳ������ʵĲ���������');
    grid on;
    
    %% 4. �ڶ�����ʹ��FFT���ƺ�ȳ�ֵ
    fprintf('\n����2: FFT���ƺ�ȳ�ֵ...\n');
    
    % ʹ������Ƶ�ʵ�������
    center_idx = round(length(wavenumber)/2);
    n_center = n1_all(center_idx);
    
    % FFT����
    thickness_init = fft_thickness_estimate(wavenumber, reflectance, n_center);
    fprintf('FFT���ƺ��: %.2f ��m\n', thickness_init);
    
    %% 5. ����������ȷ��Ϻ��
    fprintf('\n����3: ��ȷ��Ϻ��...\n');
    
    % �ĵ������ʳ�ֵ��̼�������ֵ��
    n2_init = 2.55;
    
    % �Ż���Ⱥͳĵ�������
    [thickness, n2] = fit_thickness(wavenumber, reflectance, n1_all, incident_angle, thickness_init, n2_init);
    
    fprintf('���պ��: %.2f ��m\n', thickness);
    fprintf('�ĵ�������: %.3f\n', n2);
    
    %% 6. ������Ͻ��
    plot_fitting_results(wavenumber, reflectance, n1_all, thickness, n2, incident_angle);
end

%% ������Ү�����̼���������
function n = calculate_sellmeier_n(wavelength, params)
    % Sellmeier����: n^2 = 1 + B1*��^2/(��^2-C1) + B2*��^2/(��^2-C2) + B3*��^2/(��^2-C3)
    B1 = params(1); B2 = params(2); B3 = params(3);
    C1 = params(4); C2 = params(5); C3 = params(6);
    
    lambda_sq = wavelength.^2;
    
    % ������������㣩
    term1 = B1 * lambda_sq ./ (lambda_sq - C1 + 1e-10);
    term2 = B2 * lambda_sq ./ (lambda_sq - C2 + 1e-10);  
    term3 = B3 * lambda_sq ./ (lambda_sq - C3 + 1e-10);
    
    n_squared = 1 + term1 + term2 + term3;
    
    % ȷ��nΪʵ����Ϊ��
    n = sqrt(abs(n_squared));
end

%% ���������Ү��������ʹ�ø߲�������
function params = fit_sellmeier_params(wavelength, reflectance, incident_angle)
    
    % ��ʼ�����²� [B1, B2, B3, C1, C2, C3]
    % ����̼����ĵ���ֵ
    x0 = [5.0, 0.5, 0.1, 0.01, 0.5, 20];
    
    % �����߽�
    lb = [0.1, 0.01, 0.001, 0.001, 0.01, 1];
    ub = [15, 5, 2, 1, 10, 100];
    
    % ����Ŀ�꺯������С�������ʲв�
    % �ڸ߲������򣬼�����Ӱ���С����Ҫ�������ʾ���������
    theta0 = incident_angle * pi / 180;
    
    objective = @(x) sellmeier_objective(x, wavelength, reflectance, theta0);
    
    % ʹ��lsqnonlin��ⳬ������
    options = optimoptions('lsqnonlin', ...
        'Display', 'iter', ...
        'MaxIterations', 1000, ...
        'FunctionTolerance', 1e-10, ...
        'StepTolerance', 1e-10);
    
    % ���
    [params, resnorm] = lsqnonlin(objective, x0, lb, ub, options);
    
    fprintf('��ϲв��: %.6f\n', resnorm);
end

%% ������Ү����ϵ�Ŀ�꺯��
function residuals = sellmeier_objective(params, wavelength, reflectance, theta0)
    
    % ����������
    n1 = calculate_sellmeier_n(wavelength, params);
    
    % �����ͳĵ׵�������
    n0 = 1.0;
    n2 = 2.55;  % �ĵ������ʹ���ֵ
    
    % ���������
    sin_theta1 = n0 * sin(theta0) ./ n1;
    cos_theta1 = sqrt(1 - sin_theta1.^2);
    
    % ����������ϵ�����򻯣��߲���ʱ��Ҫ���Ǳ��淴�䣩
    % sƫ��
    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    
    % pƫ��
    r01_p = (n1*cos(theta0) - n0.*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    
    % ��ƫ��ⷴ���ʣ�s��p��ƽ����
    R_surface = (abs(r01_s).^2 + abs(r01_p).^2) / 2;
    
    % �в��Ȩ�������Ͼ��ȣ�
    weight = sqrt(reflectance + 0.01);
    residuals = (R_surface - reflectance) .* weight;
end

%% FFT���ƺ��
function thickness = fft_thickness_estimate(wavenumber, reflectance, n_avg)
    
    % ȥ��ֱ������
    reflectance_ac = reflectance - mean(reflectance);
    
    % ��ֵ�����ȼ��
    k_uniform = linspace(min(wavenumber), max(wavenumber), length(wavenumber));
    r_uniform = interp1(wavenumber, reflectance_ac, k_uniform, 'spline');
    
    % FFT����
    N = length(k_uniform);
    fft_result = fft(r_uniform);
    fft_power = abs(fft_result).^2;
    
    % Ƶ���ᣨ�����
    dk = mean(diff(k_uniform));
    thickness_axis = (0:N-1) * 10000 / (2 * n_avg * N * dk);
    
    % �ҵ����壨�ų���Ƶ��
    [~, max_idx] = max(fft_power(2:round(N/2)));
    thickness = thickness_axis(max_idx + 1);
    
    % �����ں���Χ
    thickness = max(10, min(150, thickness));
end

%% ��ȷ��Ϻ��
function [thickness, n2] = fit_thickness(wavenumber, reflectance, n1, incident_angle, d_init, n2_init)
    
    % ת��Ϊ����
    theta0 = incident_angle * pi / 180;
    
    % ����Ŀ�꺯��
    objective = @(x) thickness_objective(x, wavenumber, reflectance, n1, theta0);
    
    % ��ʼֵ�ͱ߽�
    x0 = [d_init, n2_init];
    lb = [1, 2.0];
    ub = [200, 3.0];
    
    % �Ż�ѡ��
    options = optimoptions('lsqcurvefit', ...
        'Display', 'iter', ...
        'MaxIterations', 500, ...
        'FunctionTolerance', 1e-10, ...
        'StepTolerance', 1e-10);
    
    % ���
    [x_optimal, resnorm] = lsqcurvefit(@(x,xdata) model_reflectance(x, xdata, n1, theta0), ...
        x0, wavenumber, reflectance, lb, ub, options);
    
    thickness = x_optimal(1);
    n2 = x_optimal(2);
    
    fprintf('��ϲв�: %.6f\n', resnorm);
end

%% �����ϵ�Ŀ�꺯��
function residuals = thickness_objective(params, wavenumber, reflectance, n1, theta0)
    
    thickness = params(1);
    n2 = params(2);
    
    % ����ģ�ͷ�����
    R_model = model_reflectance(params, wavenumber, n1, theta0);
    
    % �в�
    residuals = R_model - reflectance;
end

%% ������ģ�ͣ��޸��汾��
function R = model_reflectance(params, wavenumber, n1, theta0)
    
    thickness = params(1);
    n2 = params(2);
    n0 = 1.0;
    
    % ȷ��n1��������
    n1 = n1(:);
    wavenumber = wavenumber(:);
    
    % ��������ǣ�������
    sin_theta1 = n0 * sin(theta0) ./ n1;
    sin_theta1 = max(-1, min(1, sin_theta1));  % ȷ������Ч��Χ
    cos_theta1 = sqrt(1 - sin_theta1.^2);
    
    % �ĵ׵�����ǣ�������
    sin_theta2 = n0 * sin(theta0) / n2;
    sin_theta2 = max(-1, min(1, sin_theta2));
    cos_theta2 = sqrt(1 - sin_theta2^2);
    
    % ������ϵ��
    % r01����������Ϊn1��������
    r01 = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    t01 = 2*n0*cos(theta0) ./ (n0*cos(theta0) + n1.*cos_theta1);
    t10 = 2*n1.*cos_theta1 ./ (n0*cos(theta0) + n1.*cos_theta1);
    
    % r12����������Ϊn1��������
    r12 = (n1.*cos_theta1 - n2*cos_theta2) ./ (n1.*cos_theta1 + n2*cos_theta2);
    
    % ��λ�������
    delta = 4 * pi * n1 .* thickness .* cos_theta1 .* wavenumber / 10000;
    
    % ˫�������湫ʽ��ȫ��ʹ�õ����㣩
    numerator = r01.^2 + r12.^2 .* t01.^2 .* t10.^2 + 2.*r01.*r12.*t01.*t10.*cos(delta);
    denominator = 1 + r01.^2 .* r12.^2 + 2.*r01.*r12.*cos(delta);
    
    R = numerator ./ denominator;
    
    % ȷ��R��������
    R = R(:);
end

%% ������Ͻ��
function plot_fitting_results(wavenumber, reflectance, n1, thickness, n2, incident_angle)
    
    theta0 = incident_angle * pi / 180;
    
    % ������ϵķ�����
    R_fitted = model_reflectance([thickness, n2], wavenumber, n1, theta0);
    
    % ����ͼ��
    figure('Position', [100, 100, 1200, 800]);
    
    % ��ͼ1�������ʶԱ�
    subplot(3,1,1);
    plot(wavenumber, reflectance*100, 'b.', 'MarkerSize', 3);
    hold on;
    plot(wavenumber, R_fitted*100, 'r-', 'LineWidth', 1.5);
    xlabel('���� (cm^{-1})');
    ylabel('������ (%)');
    title(sprintf('��������Ͻ�� (�����=%d��, ���=%.2f ��m, n_2=%.3f)', ...
        incident_angle, thickness, n2));
    legend('ʵ������', '��Ͻ��', 'Location', 'best');
    grid on;
    xlim([min(wavenumber), max(wavenumber)]);
    
    % ��ͼ2���в�ֲ�
    subplot(3,1,2);
    residuals = (reflectance - R_fitted) * 100;
    plot(wavenumber, residuals, 'g.', 'MarkerSize', 2);
    hold on;
    yline(0, 'k--', 'LineWidth', 1);
    xlabel('���� (cm^{-1})');
    ylabel('�в� (%)');
    title(sprintf('�в�ֲ� (RMSE = %.4f%%)', sqrt(mean(residuals.^2))));
    grid on;
    xlim([min(wavenumber), max(wavenumber)]);
    
    % ��ͼ3���ֲ��Ŵ�ͼ���߲�������
    subplot(3,1,3);
    high_k_mask = wavenumber > 2000;
    plot(wavenumber(high_k_mask), reflectance(high_k_mask)*100, 'b.', 'MarkerSize', 3);
    hold on;
    plot(wavenumber(high_k_mask), R_fitted(high_k_mask)*100, 'r-', 'LineWidth', 1.5);
    xlabel('���� (cm^{-1})');
    ylabel('������ (%)');
    title('�߲����������ϸ�� (>2000 cm^{-1})');
    legend('ʵ������', '��Ͻ��', 'Location', 'best');
    grid on;
    
    % ��������Ŷ�
    R_squared = 1 - sum((reflectance - R_fitted).^2) / sum((reflectance - mean(reflectance)).^2);
    
    % ����ܱ���
    sgtitle(sprintf('̼�������Ӳ������� (R? = %.4f)', R_squared));
end

%% ������������ʾ������ͳ��
function display_n_statistics(n1_all, wavenumber)
    fprintf('\n������ͳ��:\n');
    fprintf('  ��Сֵ: %.4f\n', min(n1_all));
    fprintf('  ���ֵ: %.4f\n', max(n1_all));
    fprintf('  ƽ��ֵ: %.4f\n', mean(n1_all));
    fprintf('  ��׼��: %.4f\n', std(n1_all));
    
    % �ض�������������
    specific_k = [500, 1000, 1500, 2000, 2500, 3000];
    fprintf('\n�ض�������������:\n');
    for k = specific_k
        [~, idx] = min(abs(wavenumber - k));
        if idx <= length(n1_all)
            fprintf('  %4d cm^-1: n = %.4f\n', k, n1_all(idx));
        end
    end
end