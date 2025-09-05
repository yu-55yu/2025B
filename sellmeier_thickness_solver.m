%% 主程序：塞尔迈耶尔方程拟合和外延层厚度计算
% 文件名：sellmeier_thickness_solver.m

function sellmeier_thickness_solver()
    clc; clear; close all;
    
    %% 处理附件1（10度入射角）
    fprintf('==================================================\n');
    fprintf('处理附件1数据（入射角10°）\n');
    fprintf('==================================================\n');
    
    [thickness1, n1_all_1, params1] = process_file('附件1.xlsx', 10);
    
    %% 处理附件2（15度入射角）  
    fprintf('\n==================================================\n');
    fprintf('处理附件2数据（入射角15°）\n');
    fprintf('==================================================\n');
    
    [thickness2, n1_all_2, params2] = process_file('附件2.xlsx', 15);
    
    %% 结果对比
    fprintf('\n==================================================\n');
    fprintf('结果对比分析\n');
    fprintf('==================================================\n');
    fprintf('附件1厚度: %.2f μm\n', thickness1);
    fprintf('附件2厚度: %.2f μm\n', thickness2);
    fprintf('厚度差异: %.2f μm (%.1f%%)\n', ...
        abs(thickness1-thickness2), ...
        abs(thickness1-thickness2)/mean([thickness1,thickness2])*100);
end

%% 处理单个文件
function [thickness, n1_all, sellmeier_params] = process_file(filename, incident_angle)
    
    %% 1. 读取数据
    data = readmatrix(filename);
    wavenumber = data(:, 1);  % 波数 cm^-1
    reflectance = data(:, 2) / 100;  % 反射率转换为小数
    wavelength = 10000 ./ wavenumber;  % 波长 μm
    
    %% 2. 第一步：用高波数数据拟合塞尔迈耶尔参数
    fprintf('\n步骤1: 使用高波数区域(>2000 cm^-1)拟合塞尔迈耶尔参数...\n');
    
    % 选择高波数数据
    high_k_mask = wavenumber > 2000;
    high_k_wavelength = wavelength(high_k_mask);
    high_k_reflectance = reflectance(high_k_mask);
    high_k_wavenumber = wavenumber(high_k_mask);
    
    % 拟合塞尔迈耶尔参数
    sellmeier_params = fit_sellmeier_params(high_k_wavelength, high_k_reflectance, incident_angle);
    
    fprintf('塞尔迈耶尔参数:\n');
    fprintf('  B1=%.4f, B2=%.4f, B3=%.4f\n', sellmeier_params(1), sellmeier_params(2), sellmeier_params(3));
    fprintf('  C1=%.4f, C2=%.4f, C3=%.4f\n', sellmeier_params(4), sellmeier_params(5), sellmeier_params(6));
    
    %% 3. 计算所有波长的折射率
    n1_all = calculate_sellmeier_n(wavelength, sellmeier_params);
    
    % 显示折射率统计
    display_n_statistics(n1_all, wavenumber);
    
    % 绘制折射率曲线
    figure;
    subplot(2,1,1);
    plot(wavelength, n1_all, 'b-', 'LineWidth', 1.5);
    xlabel('波长 (μm)');
    ylabel('折射率 n');
    title('外延层折射率的波长依赖性');
    grid on;
    
    subplot(2,1,2);
    plot(wavenumber, n1_all, 'r-', 'LineWidth', 1.5);
    xlabel('波数 (cm^{-1})');
    ylabel('折射率 n');
    title('外延层折射率的波数依赖性');
    grid on;
    
    %% 4. 第二步：使用FFT估计厚度初值
    fprintf('\n步骤2: FFT估计厚度初值...\n');
    
    % 使用中心频率的折射率
    center_idx = round(length(wavenumber)/2);
    n_center = n1_all(center_idx);
    
    % FFT分析
    thickness_init = fft_thickness_estimate(wavenumber, reflectance, n_center);
    fprintf('FFT估计厚度: %.2f μm\n', thickness_init);
    
    %% 5. 第三步：精确拟合厚度
    fprintf('\n步骤3: 精确拟合厚度...\n');
    
    % 衬底折射率初值（碳化硅典型值）
    n2_init = 2.55;
    
    % 优化厚度和衬底折射率
    [thickness, n2] = fit_thickness(wavenumber, reflectance, n1_all, incident_angle, thickness_init, n2_init);
    
    fprintf('最终厚度: %.2f μm\n', thickness);
    fprintf('衬底折射率: %.3f\n', n2);
    
    %% 6. 绘制拟合结果
    plot_fitting_results(wavenumber, reflectance, n1_all, thickness, n2, incident_angle);
end

%% 塞尔迈耶尔方程计算折射率
function n = calculate_sellmeier_n(wavelength, params)
    % Sellmeier方程: n^2 = 1 + B1*λ^2/(λ^2-C1) + B2*λ^2/(λ^2-C2) + B3*λ^2/(λ^2-C3)
    B1 = params(1); B2 = params(2); B3 = params(3);
    C1 = params(4); C2 = params(5); C3 = params(6);
    
    lambda_sq = wavelength.^2;
    
    % 计算各项（避免除零）
    term1 = B1 * lambda_sq ./ (lambda_sq - C1 + 1e-10);
    term2 = B2 * lambda_sq ./ (lambda_sq - C2 + 1e-10);  
    term3 = B3 * lambda_sq ./ (lambda_sq - C3 + 1e-10);
    
    n_squared = 1 + term1 + term2 + term3;
    
    % 确保n为实数且为正
    n = sqrt(abs(n_squared));
end

%% 拟合塞尔迈耶尔参数（使用高波数区域）
function params = fit_sellmeier_params(wavelength, reflectance, incident_angle)
    
    % 初始参数猜测 [B1, B2, B3, C1, C2, C3]
    % 基于碳化硅的典型值
    x0 = [5.0, 0.5, 0.1, 0.01, 0.5, 20];
    
    % 参数边界
    lb = [0.1, 0.01, 0.001, 0.001, 0.01, 1];
    ub = [15, 5, 2, 1, 10, 100];
    
    % 定义目标函数：最小化反射率残差
    % 在高波数区域，假设厚度影响较小，主要由折射率决定反射率
    theta0 = incident_angle * pi / 180;
    
    objective = @(x) sellmeier_objective(x, wavelength, reflectance, theta0);
    
    % 使用lsqnonlin求解超定方程
    options = optimoptions('lsqnonlin', ...
        'Display', 'iter', ...
        'MaxIterations', 1000, ...
        'FunctionTolerance', 1e-10, ...
        'StepTolerance', 1e-10);
    
    % 求解
    [params, resnorm] = lsqnonlin(objective, x0, lb, ub, options);
    
    fprintf('拟合残差范数: %.6f\n', resnorm);
end

%% 塞尔迈耶尔拟合的目标函数
function residuals = sellmeier_objective(params, wavelength, reflectance, theta0)
    
    % 计算折射率
    n1 = calculate_sellmeier_n(wavelength, params);
    
    % 空气和衬底的折射率
    n0 = 1.0;
    n2 = 2.55;  % 衬底折射率估计值
    
    % 计算折射角
    sin_theta1 = n0 * sin(theta0) ./ n1;
    cos_theta1 = sqrt(1 - sin_theta1.^2);
    
    % 菲涅尔反射系数（简化：高波数时主要考虑表面反射）
    % s偏振
    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    
    % p偏振
    r01_p = (n1*cos(theta0) - n0.*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    
    % 非偏振光反射率（s和p的平均）
    R_surface = (abs(r01_s).^2 + abs(r01_p).^2) / 2;
    
    % 残差（加权以提高拟合精度）
    weight = sqrt(reflectance + 0.01);
    residuals = (R_surface - reflectance) .* weight;
end

%% FFT估计厚度
function thickness = fft_thickness_estimate(wavenumber, reflectance, n_avg)
    
    % 去除直流分量
    reflectance_ac = reflectance - mean(reflectance);
    
    % 插值到均匀间隔
    k_uniform = linspace(min(wavenumber), max(wavenumber), length(wavenumber));
    r_uniform = interp1(wavenumber, reflectance_ac, k_uniform, 'spline');
    
    % FFT分析
    N = length(k_uniform);
    fft_result = fft(r_uniform);
    fft_power = abs(fft_result).^2;
    
    % 频率轴（厚度域）
    dk = mean(diff(k_uniform));
    thickness_axis = (0:N-1) * 10000 / (2 * n_avg * N * dk);
    
    % 找到主峰（排除零频）
    [~, max_idx] = max(fft_power(2:round(N/2)));
    thickness = thickness_axis(max_idx + 1);
    
    % 限制在合理范围
    thickness = max(10, min(150, thickness));
end

%% 精确拟合厚度
function [thickness, n2] = fit_thickness(wavenumber, reflectance, n1, incident_angle, d_init, n2_init)
    
    % 转换为弧度
    theta0 = incident_angle * pi / 180;
    
    % 定义目标函数
    objective = @(x) thickness_objective(x, wavenumber, reflectance, n1, theta0);
    
    % 初始值和边界
    x0 = [d_init, n2_init];
    lb = [1, 2.0];
    ub = [200, 3.0];
    
    % 优化选项
    options = optimoptions('lsqcurvefit', ...
        'Display', 'iter', ...
        'MaxIterations', 500, ...
        'FunctionTolerance', 1e-10, ...
        'StepTolerance', 1e-10);
    
    % 拟合
    [x_optimal, resnorm] = lsqcurvefit(@(x,xdata) model_reflectance(x, xdata, n1, theta0), ...
        x0, wavenumber, reflectance, lb, ub, options);
    
    thickness = x_optimal(1);
    n2 = x_optimal(2);
    
    fprintf('拟合残差: %.6f\n', resnorm);
end

%% 厚度拟合的目标函数
function residuals = thickness_objective(params, wavenumber, reflectance, n1, theta0)
    
    thickness = params(1);
    n2 = params(2);
    
    % 计算模型反射率
    R_model = model_reflectance(params, wavenumber, n1, theta0);
    
    % 残差
    residuals = R_model - reflectance;
end

%% 反射率模型（修复版本）
function R = model_reflectance(params, wavenumber, n1, theta0)
    
    thickness = params(1);
    n2 = params(2);
    n0 = 1.0;
    
    % 确保n1是列向量
    n1 = n1(:);
    wavenumber = wavenumber(:);
    
    % 计算折射角（向量）
    sin_theta1 = n0 * sin(theta0) ./ n1;
    sin_theta1 = max(-1, min(1, sin_theta1));  % 确保在有效范围
    cos_theta1 = sqrt(1 - sin_theta1.^2);
    
    % 衬底的折射角（标量）
    sin_theta2 = n0 * sin(theta0) / n2;
    sin_theta2 = max(-1, min(1, sin_theta2));
    cos_theta2 = sqrt(1 - sin_theta2^2);
    
    % 菲涅尔系数
    % r01是向量（因为n1是向量）
    r01 = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    t01 = 2*n0*cos(theta0) ./ (n0*cos(theta0) + n1.*cos_theta1);
    t10 = 2*n1.*cos_theta1 ./ (n0*cos(theta0) + n1.*cos_theta1);
    
    % r12是向量（因为n1是向量）
    r12 = (n1.*cos_theta1 - n2*cos_theta2) ./ (n1.*cos_theta1 + n2*cos_theta2);
    
    % 相位差（向量）
    delta = 4 * pi * n1 .* thickness .* cos_theta1 .* wavenumber / 10000;
    
    % 双光束干涉公式（全部使用点运算）
    numerator = r01.^2 + r12.^2 .* t01.^2 .* t10.^2 + 2.*r01.*r12.*t01.*t10.*cos(delta);
    denominator = 1 + r01.^2 .* r12.^2 + 2.*r01.*r12.*cos(delta);
    
    R = numerator ./ denominator;
    
    % 确保R是列向量
    R = R(:);
end

%% 绘制拟合结果
function plot_fitting_results(wavenumber, reflectance, n1, thickness, n2, incident_angle)
    
    theta0 = incident_angle * pi / 180;
    
    % 计算拟合的反射率
    R_fitted = model_reflectance([thickness, n2], wavenumber, n1, theta0);
    
    % 创建图形
    figure('Position', [100, 100, 1200, 800]);
    
    % 子图1：反射率对比
    subplot(3,1,1);
    plot(wavenumber, reflectance*100, 'b.', 'MarkerSize', 3);
    hold on;
    plot(wavenumber, R_fitted*100, 'r-', 'LineWidth', 1.5);
    xlabel('波数 (cm^{-1})');
    ylabel('反射率 (%)');
    title(sprintf('反射率拟合结果 (入射角=%d°, 厚度=%.2f μm, n_2=%.3f)', ...
        incident_angle, thickness, n2));
    legend('实测数据', '拟合结果', 'Location', 'best');
    grid on;
    xlim([min(wavenumber), max(wavenumber)]);
    
    % 子图2：残差分布
    subplot(3,1,2);
    residuals = (reflectance - R_fitted) * 100;
    plot(wavenumber, residuals, 'g.', 'MarkerSize', 2);
    hold on;
    yline(0, 'k--', 'LineWidth', 1);
    xlabel('波数 (cm^{-1})');
    ylabel('残差 (%)');
    title(sprintf('残差分布 (RMSE = %.4f%%)', sqrt(mean(residuals.^2))));
    grid on;
    xlim([min(wavenumber), max(wavenumber)]);
    
    % 子图3：局部放大图（高波数区域）
    subplot(3,1,3);
    high_k_mask = wavenumber > 2000;
    plot(wavenumber(high_k_mask), reflectance(high_k_mask)*100, 'b.', 'MarkerSize', 3);
    hold on;
    plot(wavenumber(high_k_mask), R_fitted(high_k_mask)*100, 'r-', 'LineWidth', 1.5);
    xlabel('波数 (cm^{-1})');
    ylabel('反射率 (%)');
    title('高波数区域拟合细节 (>2000 cm^{-1})');
    legend('实测数据', '拟合结果', 'Location', 'best');
    grid on;
    
    % 计算拟合优度
    R_squared = 1 - sum((reflectance - R_fitted).^2) / sum((reflectance - mean(reflectance)).^2);
    
    % 添加总标题
    sgtitle(sprintf('碳化硅外延层测量结果 (R? = %.4f)', R_squared));
end

%% 辅助函数：显示折射率统计
function display_n_statistics(n1_all, wavenumber)
    fprintf('\n折射率统计:\n');
    fprintf('  最小值: %.4f\n', min(n1_all));
    fprintf('  最大值: %.4f\n', max(n1_all));
    fprintf('  平均值: %.4f\n', mean(n1_all));
    fprintf('  标准差: %.4f\n', std(n1_all));
    
    % 特定波数的折射率
    specific_k = [500, 1000, 1500, 2000, 2500, 3000];
    fprintf('\n特定波数的折射率:\n');
    for k = specific_k
        [~, idx] = min(abs(wavenumber - k));
        if idx <= length(n1_all)
            fprintf('  %4d cm^-1: n = %.4f\n', k, n1_all(idx));
        end
    end
end