%% 主函数：碳化硅外延层厚度测量
% 文件名：main_sic_thickness.m
function main_sic_thickness()
    clc; clear; close all;
    
    %% 处理附件1（10度入射角）
    fprintf('==================================================\n');
    fprintf('处理附件1数据（入射角10°）\n');
    fprintf('==================================================\n');
    
    [thickness1, n2_1, params1, rmse1] = process_data('附件1.xlsx', 10);
    
    fprintf('\n附件1结果:\n');
    fprintf('外延层厚度: %.2f μm\n', thickness1);
    fprintf('衬底折射率: %.3f\n', n2_1);
    fprintf('塞尔迈耶尔参数:\n');
    fprintf('  B1=%.3f, B2=%.3f, B3=%.3f\n', params1(1), params1(2), params1(3));
    fprintf('  C1=%.3f, C2=%.3f, C3=%.3f\n', params1(4), params1(5), params1(6));
    fprintf('RMSE: %.4f\n', rmse1);
    
    %% 处理附件2（15度入射角）
    fprintf('\n==================================================\n');
    fprintf('处理附件2数据（入射角15°）\n');
    fprintf('==================================================\n');
    
    [thickness2, n2_2, params2, rmse2] = process_data('附件2.xlsx', 15);
    
    fprintf('\n附件2结果:\n');
    fprintf('外延层厚度: %.2f μm\n', thickness2);
    fprintf('衬底折射率: %.3f\n', n2_2);
    fprintf('塞尔迈耶尔参数:\n');
    fprintf('  B1=%.3f, B2=%.3f, B3=%.3f\n', params2(1), params2(2), params2(3));
    fprintf('  C1=%.3f, C2=%.3f, C3=%.3f\n', params2(4), params2(5), params2(6));
    fprintf('RMSE: %.4f\n', rmse2);
    
    %% 可靠性分析
    fprintf('\n==================================================\n');
    fprintf('可靠性分析:\n');
    fprintf('==================================================\n');
    
    thickness_diff = abs(thickness1 - thickness2);
    relative_diff = thickness_diff / ((thickness1 + thickness2) / 2) * 100;
    
    fprintf('两种角度测得的厚度差: %.2f μm\n', thickness_diff);
    fprintf('相对偏差: %.1f%%\n', relative_diff);
    
    if thickness_diff < 1.0
        fprintf('结果可靠性: 优秀（厚度差<1μm）\n');
    elseif thickness_diff < 2.0
        fprintf('结果可靠性: 良好（厚度差<2μm）\n');
    else
        fprintf('结果可靠性: 需要进一步验证\n');
    end
end

%% 处理单个数据文件
function [thickness, n2, sellmeier_params, rmse] = process_data(filename, incident_angle)
    % 读取数据
    data = readmatrix(filename);
    wavenumber = data(:, 1);  % 波数 cm^-1
    reflectance = data(:, 2) / 100;  % 反射率转换为小数
    
    % 波长(μm) = 10000/波数
    wavelength = 10000 ./ wavenumber;
    
    % 入射角转换为弧度
    theta0 = incident_angle * pi / 180;
    
    % 使用FFT估计厚度初值
    thickness_init = estimate_thickness_fft(wavenumber, reflectance);
    fprintf('FFT估计厚度初值: %.2f μm\n', thickness_init);
    
    % 设置优化参数
    % 参数顺序: [B1, B2, B3, C1, C2, C3, thickness, n2]
    x0 = [2.0, 0.5, 0.1, 0.01, 0.1, 10, thickness_init, 2.55];
    lb = [0.1, 0.01, 0.001, 0.001, 0.01, 1, 1, 2.0];
    ub = [10, 5, 2, 1, 10, 100, 200, 3.0];
    
    % 定义目标函数
    objective = @(x) sum(calculate_residuals(x, wavelength, wavenumber, reflectance, theta0).^2);
    
    % 全局优化选项
    fprintf('开始全局优化...\n');
    
    % 使用遗传算法进行全局优化
    options_ga = optimoptions('ga', ...
        'PopulationSize', 100, ...
        'MaxGenerations', 50, ...
        'Display', 'iter', ...
        'UseParallel', true, ...
        'FunctionTolerance', 1e-8);
    
    [x_ga, ~] = ga(objective, 8, [], [], [], [], lb, ub, [], options_ga);
    
    % 使用局部优化进一步精细化
    fprintf('\n开始局部优化...\n');
    options_local = optimoptions('fmincon', ...
        'Display', 'iter', ...
        'Algorithm', 'interior-point', ...
        'MaxIterations', 500, ...
        'OptimalityTolerance', 1e-10, ...
        'StepTolerance', 1e-10);
    
    [x_optimal, fval] = fmincon(objective, x_ga, [], [], [], [], lb, ub, [], options_local);
    
    % 提取结果
    sellmeier_params = x_optimal(1:6);
    thickness = x_optimal(7);
    n2 = x_optimal(8);
    
    % 计算RMSE
    residuals = calculate_residuals(x_optimal, wavelength, wavenumber, reflectance, theta0);
    rmse = sqrt(mean(residuals.^2));
    
    % 绘制结果
    plot_results(wavenumber, reflectance, x_optimal, wavelength, theta0, incident_angle);
end

%% FFT估计厚度
function thickness_est = estimate_thickness_fft(wavenumber, reflectance)
    % 中心化数据
    reflectance_centered = reflectance - mean(reflectance);
    
    % FFT分析
    N = length(wavenumber);
    fft_result = fft(reflectance_centered);
    
    % 计算频率
    delta_k = mean(diff(wavenumber));
    freq = (0:N-1) / (N * delta_k);
    
    % 找到主频率（排除DC分量）
    fft_power = abs(fft_result).^2;
    [~, max_idx] = max(fft_power(2:floor(N/2)));
    main_freq = freq(max_idx + 1);
    
    % 估计厚度 d ≈ 10000/(2n*主频率)
    n_avg = 2.6;  % 碳化硅典型折射率
    thickness_est = 10000 / (2 * n_avg * main_freq);
    
    % 限制在合理范围内
    thickness_est = max(10, min(150, thickness_est));
end

%% 计算塞尔迈耶尔折射率
function n = sellmeier(wavelength, params)
    B1 = params(1); B2 = params(2); B3 = params(3);
    C1 = params(4); C2 = params(5); C3 = params(6);
    
    lambda_sq = wavelength.^2;
    
    % 避免除零
    term1 = B1 * lambda_sq ./ (lambda_sq - C1 + 1e-10);
    term2 = B2 * lambda_sq ./ (lambda_sq - C2 + 1e-10);
    term3 = B3 * lambda_sq ./ (lambda_sq - C3 + 1e-10);
    
    n_squared = 1 + term1 + term2 + term3;
    n = sqrt(abs(n_squared));
end

%% 计算反射率模型
function R = reflectance_model(params, wavelength, wavenumber, theta0)
    % 提取参数
    sellmeier_params = params(1:6);
    thickness = params(7);
    n2 = params(8);
    
    % 计算外延层折射率
    n1 = sellmeier(wavelength, sellmeier_params);
    n0 = 1.0;  % 空气折射率
    
    % 计算折射角（斯涅尔定律）
    sin_theta1 = n0 * sin(theta0) ./ n1;
    sin_theta2 = n0 * sin(theta0) / n2;
    
    % 确保值在有效范围内
    sin_theta1 = max(-1, min(1, sin_theta1));
    sin_theta2 = max(-1, min(1, sin_theta2));
    
    cos_theta0 = cos(theta0);
    cos_theta1 = sqrt(1 - sin_theta1.^2);
    cos_theta2 = sqrt(1 - sin_theta2^2);
    
    % 计算菲涅尔系数（s偏振）
    r01 = (n0*cos_theta0 - n1.*cos_theta1) ./ (n0*cos_theta0 + n1.*cos_theta1);
    t01 = 2*n0*cos_theta0 ./ (n0*cos_theta0 + n1.*cos_theta1);
    t10 = 2*n1.*cos_theta1 ./ (n0*cos_theta0 + n1.*cos_theta1);
    r12 = (n1.*cos_theta1 - n2*cos_theta2) ./ (n1.*cos_theta1 + n2*cos_theta2);
    
    % 相位差 δ = 4πn₁d·cos(θ₁)·ν/10000
    delta = 4 * pi * n1 .* thickness .* cos_theta1 .* wavenumber / 10000;
    
    % 计算反射率
    numerator = r01.^2 + r12^2 * t01.^2 .* t10.^2 + 2*r01.*r12.*t01.*t10.*cos(delta);
    denominator = 1 + r01.^2 * r12^2 + 2*r01.*r12.*cos(delta);
    
    R = numerator ./ denominator;
end

%% 计算残差
function residuals = calculate_residuals(params, wavelength, wavenumber, reflectance, theta0)
    R_model = reflectance_model(params, wavelength, wavenumber, theta0);
    % 加权残差
    weights = sqrt(reflectance + 0.01);
    residuals = (R_model - reflectance) .* weights;
end

%% 绘制结果
function plot_results(wavenumber, reflectance, params, wavelength, theta0, incident_angle)
    % 计算拟合的反射率
    R_fitted = reflectance_model(params, wavelength, wavenumber, theta0);
    
    % 提取参数
    thickness = params(7);
    n2 = params(8);
    
    % 创建图形
    figure('Position', [100, 100, 1000, 600]);
    
    % 子图1：反射率对比
    subplot(2, 1, 1);
    plot(wavenumber, reflectance, 'b.', 'MarkerSize', 2);
    hold on;
    plot(wavenumber, R_fitted, 'r-', 'LineWidth', 1.5);
    xlabel('波数 (cm^{-1})');
    ylabel('反射率');
    title(sprintf('反射率拟合结果 (入射角=%d°, 厚度=%.2f μm, n_2=%.3f)', ...
        incident_angle, thickness, n2));
    legend('实测数据', '拟合结果', 'Location', 'best');
    grid on;
    grid minor;
    
    % 子图2：残差分布
    subplot(2, 1, 2);
    residuals = reflectance - R_fitted;
    plot(wavenumber, residuals, 'g.', 'MarkerSize', 2);
    hold on;
    yline(0, 'k--', 'LineWidth', 1);
    xlabel('波数 (cm^{-1})');
    ylabel('残差');
    title(sprintf('残差分布 (RMSE = %.4f)', sqrt(mean(residuals.^2))));
    grid on;
    grid minor;
    
    % 调整布局
    sgtitle(sprintf('碳化硅外延层厚度测量结果 (入射角 %d°)', incident_angle));
end

%% 获取所有波长对应的折射率
function n1_all = get_n1_all(wavelength, sellmeier_params)
    n1_all = sellmeier(wavelength, sellmeier_params);
end

%% 额外的辅助函数：批量处理和分析
function batch_analysis()
    % 此函数用于批量处理多个文件或进行参数敏感性分析
    
    % 角度扫描分析
    angles = 5:5:30;
    thicknesses = zeros(length(angles), 1);
    
    for i = 1:length(angles)
        % 这里假设有对应角度的数据文件
        filename = sprintf('data_angle_%d.xlsx', angles(i));
        if exist(filename, 'file')
            [thickness, ~, ~, ~] = process_data(filename, angles(i));
            thicknesses(i) = thickness;
        end
    end
    
    % 绘制角度依赖性
    figure;
    plot(angles, thicknesses, 'bo-', 'LineWidth', 2);
    xlabel('入射角 (度)');
    ylabel('测量厚度 (μm)');
    title('厚度测量的角度依赖性');
    grid on;
end