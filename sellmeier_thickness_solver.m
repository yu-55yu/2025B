%% 主程序：塞尔迈耶尔方程拟合和外延层厚度计算
% 文件名：sellmeier_thickness_solver_improved.m

function sellmeier_thickness_solver_improved()
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
    
    %% 结果对比与验证
    fprintf('\n==================================================\n');
    fprintf('结果对比分析\n');
    fprintf('==================================================\n');
    fprintf('附件1厚度: %.2f μm\n', thickness1);
    fprintf('附件2厚度: %.2f μm\n', thickness2);
    fprintf('厚度差异: %.2f μm (%.1f%%)\n', ...
        abs(thickness1-thickness2), ...
        abs(thickness1-thickness2)/mean([thickness1,thickness2])*100);
    
    % 结果可靠性分析
    analyze_reliability(thickness1, thickness2, params1, params2);
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
    
    % 选择高波数数据（只用于拟合参数）
    high_k_mask = wavenumber > 2000;
    high_k_wavelength = wavelength(high_k_mask);
    high_k_reflectance = reflectance(high_k_mask);
    
    % 拟合塞尔迈耶尔参数（改进版）
    sellmeier_params = fit_sellmeier_params_improved(high_k_wavelength, high_k_reflectance, incident_angle);
    
    fprintf('塞尔迈耶尔参数:\n');
    fprintf('  B1=%.4f, B2=%.4f, B3=%.4f\n', sellmeier_params(1), sellmeier_params(2), sellmeier_params(3));
    fprintf('  C1=%.4f, C2=%.4f, C3=%.4f\n', sellmeier_params(4), sellmeier_params(5), sellmeier_params(6));
    
    %% 3. 计算所有波长的折射率
    n1_all = calculate_sellmeier_n(wavelength, sellmeier_params);
    
    % 显示折射率统计
    display_n_statistics(n1_all, wavenumber);
    
    % 绘制折射率曲线（全波段）
    plot_refractive_index(wavelength, wavenumber, n1_all, high_k_mask);
    
    %% 4. 第二步：使用FFT估计厚度初值
    fprintf('\n步骤2: FFT估计厚度初值...\n');
    
    % 使用平均折射率
    n_avg = mean(n1_all);
    
    % FFT分析（改进版）
    thickness_init = fft_thickness_estimate_improved(wavenumber, reflectance, n_avg);
    fprintf('FFT估计厚度: %.2f μm\n', thickness_init);
    
    %% 5. 第三步：精确拟合厚度
    fprintf('\n步骤3: 精确拟合厚度...\n');
    
    % 衬底折射率初值（碳化硅典型值）
    n2_init = 2.55;
    
    % 优化厚度和衬底折射率（改进版）
    [thickness, n2, fit_quality] = fit_thickness_improved(wavenumber, reflectance, n1_all, incident_angle, thickness_init, n2_init);
    
    fprintf('最终厚度: %.2f μm\n', thickness);
    fprintf('衬底折射率: %.3f\n', n2);
    fprintf('拟合优度 R?: %.4f\n', fit_quality);
    
    %% 6. 绘制拟合结果
    plot_fitting_results_improved(wavenumber, reflectance, n1_all, thickness, n2, incident_angle);
end

%% 改进的塞尔迈耶尔参数拟合
function params = fit_sellmeier_params_improved(wavelength, reflectance, incident_angle)
    
    % 多组初始值尝试，选择最优结果
    initial_guesses = [
        [5.0, 0.5, 0.1, 0.01, 0.5, 20];    % 基础猜测
        [6.5, 0.3, 0.05, 0.02, 0.3, 15];   % 碳化硅典型值
        [4.0, 0.8, 0.2, 0.005, 0.8, 25];   % 替代猜测
    ];
    
    % 参数边界（更宽松）
    lb = [0.1, 0.001, 0.0001, 0.0001, 0.001, 0.1];
    ub = [20, 10, 5, 2, 20, 200];
    
    theta0 = incident_angle * pi / 180;
    
    best_params = [];
    best_resnorm = inf;
    
    % 尝试不同初始值
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
    fprintf('最优拟合残差: %.6f\n', best_resnorm);
end

%% 改进的塞尔迈耶尔目标函数
function residuals = sellmeier_objective_improved(params, wavelength, reflectance, theta0)
    
    % 计算折射率
    n1 = calculate_sellmeier_n(wavelength, params);
    
    % 确保折射率在合理范围
    if any(n1 < 1.5) || any(n1 > 4.0) || any(imag(n1) ~= 0)
        residuals = ones(size(reflectance)) * 1e6;  % 惩罚项
        return;
    end
    
    % 空气和衬底的折射率
    n0 = 1.0;
    
    % 计算折射角（斯涅尔定律）
    sin_theta1 = n0 * sin(theta0) ./ n1;
    cos_theta1 = real(sqrt(1 - sin_theta1.^2));
    
    % 非偏振光的菲涅尔反射系数
    % s偏振
    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    
    % p偏振
    r01_p = (n1*cos(theta0) - n0.*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    
    % 非偏振光反射率
    R_surface = (abs(r01_s).^2 + abs(r01_p).^2) / 2;
    
    % 残差（使用相对误差加权）
    weight = 1 ./ (reflectance + 0.01);  % 避免除零
    residuals = (R_surface - reflectance) .* sqrt(weight);
end

%% 改进的FFT厚度估计
function thickness = fft_thickness_estimate_improved(wavenumber, reflectance, n_avg)
    
    % 数据预处理
    reflectance_smooth = smooth(reflectance, 5);  % 平滑处理
    reflectance_ac = reflectance_smooth - mean(reflectance_smooth);
    
    % 插值到更密的均匀网格
    N = 2^nextpow2(4*length(wavenumber));  % 使用2的幂次以提高FFT效率
    k_uniform = linspace(min(wavenumber), max(wavenumber), N);
    r_uniform = interp1(wavenumber, reflectance_ac, k_uniform, 'pchip');
    
    % 加窗减少频谱泄露
    window = hann(N);
    r_windowed = r_uniform(:) .* window(:);
    
    % FFT分析
    fft_result = fft(r_windowed);
    fft_power = abs(fft_result).^2;
    
    % 频率轴（厚度域）
    dk = mean(diff(k_uniform));
    thickness_axis = (0:N-1) * 10000 / (2 * n_avg * N * dk);
    
    % 找到主峰（排除零频和噪声）
    search_range = find(thickness_axis > 5 & thickness_axis < 200);
    [~, max_idx] = max(fft_power(search_range));
    thickness = thickness_axis(search_range(max_idx));
end

%% 改进的厚度拟合
function [thickness, n2, R_squared] = fit_thickness_improved(wavenumber, reflectance, n1, incident_angle, d_init, n2_init)
    
    theta0 = incident_angle * pi / 180;
    
    % 定义目标函数
    model_func = @(x, xdata) model_reflectance_vectorized(x, xdata, n1, theta0);
    
    % 初始值和边界
    x0 = [d_init, n2_init];
    lb = [0.5, 2.0];   % 更宽的范围
    ub = [300, 3.5];
    
    % 优化选项
    options = optimoptions('lsqcurvefit', ...
        'Display', 'iter', ...
        'MaxIterations', 1000, ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'Algorithm', 'trust-region-reflective', ...
        'UseParallel', true);  % 并行计算
    
    % 拟合
    [x_optimal, resnorm] = lsqcurvefit(model_func, x0, wavenumber, reflectance, lb, ub, options);
    
    thickness = x_optimal(1);
    n2 = x_optimal(2);
    
    % 计算拟合优度
    R_fitted = model_func(x_optimal, wavenumber);
    SS_tot = sum((reflectance - mean(reflectance)).^2);
    SS_res = sum((reflectance - R_fitted).^2);
    R_squared = 1 - SS_res/SS_tot;
    
    fprintf('最终拟合残差: %.6f\n', resnorm);
end

%% 向量化的反射率模型
function R = model_reflectance_vectorized(params, wavenumber, n1, theta0)
    
    thickness = params(1);
    n2 = params(2);
    n0 = 1.0;
    
    % 确保输入是列向量
    n1 = n1(:);
    wavenumber = wavenumber(:);
    
    % 计算各层折射角
    sin_theta1 = n0 * sin(theta0) ./ n1;
    cos_theta1 = real(sqrt(1 - sin_theta1.^2));
    
    sin_theta2 = n0 * sin(theta0) / n2;
    cos_theta2 = real(sqrt(1 - sin_theta2^2));
    
    % 非偏振光的综合反射系数
    % s偏振
    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    r12_s = (n1.*cos_theta1 - n2*cos_theta2) ./ (n1.*cos_theta1 + n2*cos_theta2);
    
    % p偏振
    r01_p = (n1*cos(theta0) - n0.*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    r12_p = (n2*cos_theta1 - n1.*cos_theta2) ./ (n2*cos_theta1 + n1.*cos_theta2);
    
    % 相位差
    delta = 4 * pi * n1 .* thickness .* cos_theta1 .* wavenumber / 10000;
    
    % 双光束干涉公式（s和p的平均）
    R_s = abs((r01_s + r12_s .* exp(1i*delta)) ./ (1 + r01_s .* r12_s .* exp(1i*delta))).^2;
    R_p = abs((r01_p + r12_p .* exp(1i*delta)) ./ (1 + r01_p .* r12_p .* exp(1i*delta))).^2;
    
    R = (R_s + R_p) / 2;
    R = real(R(:));  % 确保输出是实数列向量
end

%% 塞尔迈耶尔方程计算折射率
function n = calculate_sellmeier_n(wavelength, params)
    % Sellmeier方程: n^2 = 1 + Σ[Bi*λ^2/(λ^2-Ci)]
    B1 = params(1); B2 = params(2); B3 = params(3);
    C1 = params(4); C2 = params(5); C3 = params(6);
    
    lambda_sq = wavelength.^2;
    
    % 计算各项（改进数值稳定性）
    eps = 1e-10;  % 避免除零
    term1 = B1 * lambda_sq ./ (lambda_sq - C1 + eps);
    term2 = B2 * lambda_sq ./ (lambda_sq - C2 + eps);  
    term3 = B3 * lambda_sq ./ (lambda_sq - C3 + eps);
    
    % 确保每项都是正值贡献
    term1 = max(0, real(term1));
    term2 = max(0, real(term2));
    term3 = max(0, real(term3));
    
    n_squared = 1 + term1 + term2 + term3;
    
    % 确保n为实数且合理
    n = real(sqrt(n_squared));
    n = max(1.0, min(4.0, n));  % 限制在合理范围
end

%% 绘制折射率曲线
function plot_refractive_index(wavelength, wavenumber, n1_all, high_k_mask)
    
    figure('Position', [100, 100, 1200, 600]);
    
    % 子图1：折射率vs波长（全波段）
    subplot(1,2,1);
    plot(wavelength, n1_all, 'b-', 'LineWidth', 2);
    hold on;
    plot(wavelength(high_k_mask), n1_all(high_k_mask), 'r.', 'MarkerSize', 8);
    xlabel('波长 (μm)');
    ylabel('折射率 n');
    title('外延层折射率色散曲线');
    legend('全波段折射率', '高波数拟合区域', 'Location', 'best');
    grid on;
    ylim([min(n1_all)*0.98, max(n1_all)*1.02]);
    
    % 子图2：折射率vs波数（全波段）
    subplot(1,2,2);
    plot(wavenumber, n1_all, 'b-', 'LineWidth', 2);
    hold on;
    plot(wavenumber(high_k_mask), n1_all(high_k_mask), 'r.', 'MarkerSize', 8);
    xline(2000, 'k--', '拟合边界', 'LineWidth', 1.5);
    xlabel('波数 (cm^{-1})');
    ylabel('折射率 n');
    title('外延层折射率（波数域）');
    legend('全波段折射率', '高波数拟合区域', 'Location', 'best');
    grid on;
    xlim([min(wavenumber), max(wavenumber)]);
    
    sgtitle('塞尔迈耶尔方程拟合结果');
end

%% 改进的拟合结果绘图
function plot_fitting_results_improved(wavenumber, reflectance, n1, thickness, n2, incident_angle)
    
    theta0 = incident_angle * pi / 180;
    
    % 计算拟合的反射率
    R_fitted = model_reflectance_vectorized([thickness, n2], wavenumber, n1, theta0);
    
    % 创建综合图形
    figure('Position', [100, 100, 1400, 900]);
    
    % 子图1：全波段反射率对比
    subplot(3,2,[1,2]);
    plot(wavenumber, reflectance*100, 'b.', 'MarkerSize', 2);
    hold on;
    plot(wavenumber, R_fitted*100, 'r-', 'LineWidth', 1.5);
    xlabel('波数 (cm^{-1})');
    ylabel('反射率 (%)');
    title(sprintf('反射率拟合结果 (θ=%d°, d=%.2f μm, n?=%.3f)', ...
        incident_angle, thickness, n2));
    legend('实测数据', '拟合结果', 'Location', 'best');
    grid on;
    xlim([min(wavenumber), max(wavenumber)]);
    
    % 子图2：残差分析
    subplot(3,2,3);
    residuals = (reflectance - R_fitted) * 100;
    plot(wavenumber, residuals, 'g.', 'MarkerSize', 2);
    hold on;
    yline(0, 'k--', 'LineWidth', 1);
    yline(std(residuals), 'r--', '+σ');
    yline(-std(residuals), 'r--', '-σ');
    xlabel('波数 (cm^{-1})');
    ylabel('残差 (%)');
    title(sprintf('残差分布 (RMSE = %.4f%%)', sqrt(mean(residuals.^2))));
    grid on;
    xlim([min(wavenumber), max(wavenumber)]);
    
    % 子图3：残差直方图
    subplot(3,2,4);
    histogram(residuals, 30, 'Normalization', 'probability');
    xlabel('残差 (%)');
    ylabel('概率');
    title('残差分布直方图');
    grid on;
    
    % 子图4：高波数区域细节
    subplot(3,2,5);
    high_k_mask = wavenumber > 2000;
    plot(wavenumber(high_k_mask), reflectance(high_k_mask)*100, 'b.', 'MarkerSize', 3);
    hold on;
    plot(wavenumber(high_k_mask), R_fitted(high_k_mask)*100, 'r-', 'LineWidth', 1.5);
    xlabel('波数 (cm^{-1})');
    ylabel('反射率 (%)');
    title('高波数区域拟合细节 (>2000 cm^{-1})');
    legend('实测数据', '拟合结果', 'Location', 'best');
    grid on;
    
    % 子图5：低波数区域细节
    subplot(3,2,6);
    low_k_mask = wavenumber < 1000;
    plot(wavenumber(low_k_mask), reflectance(low_k_mask)*100, 'b.', 'MarkerSize', 3);
    hold on;
    plot(wavenumber(low_k_mask), R_fitted(low_k_mask)*100, 'r-', 'LineWidth', 1.5);
    xlabel('波数 (cm^{-1})');
    ylabel('反射率 (%)');
    title('低波数区域拟合细节 (<1000 cm^{-1})');
    legend('实测数据', '拟合结果', 'Location', 'best');
    grid on;
    
    % 计算并显示拟合优度指标
    R_squared = 1 - sum((reflectance - R_fitted).^2) / sum((reflectance - mean(reflectance)).^2);
    MAE = mean(abs(residuals));
    
    % 添加总标题
    sgtitle(sprintf('碳化硅外延层厚度测量结果 (R? = %.4f, MAE = %.4f%%)', R_squared, MAE));
end

%% 显示折射率统计
function display_n_statistics(n1_all, wavenumber)
    fprintf('\n折射率统计信息:\n');
    fprintf('  范围: [%.4f, %.4f]\n', min(n1_all), max(n1_all));
    fprintf('  平均值: %.4f\n', mean(n1_all));
    fprintf('  标准差: %.4f\n', std(n1_all));
    fprintf('  中位数: %.4f\n', median(n1_all));
    
    % 特定波数点的折射率
    specific_k = [500, 1000, 1500, 2000, 2500, 3000, 3500];
    fprintf('\n特定波数点的折射率:\n');
    for k = specific_k
        [~, idx] = min(abs(wavenumber - k));
        if idx <= length(n1_all)
            fprintf('  %4d cm^-1: n = %.4f (λ = %.2f μm)\n', ...
                k, n1_all(idx), 10000/k);
        end
    end
end

%% 结果可靠性分析
function analyze_reliability(thickness1, thickness2, params1, params2)
    
    fprintf('\n=== 结果可靠性分析 ===\n');
    
    % 1. 厚度一致性检验
    thickness_diff_percent = abs(thickness1-thickness2)/mean([thickness1,thickness2])*100;
    if thickness_diff_percent < 5
        fprintf('? 厚度一致性良好 (差异 < 5%%)\n');
    else
        fprintf('? 厚度差异较大 (%.1f%%)，需要进一步检查\n', thickness_diff_percent);
    end
    
    % 2. 塞尔迈耶尔参数合理性
    fprintf('\n塞尔迈耶尔参数对比:\n');
    param_names = {'B1', 'B2', 'B3', 'C1', 'C2', 'C3'};
    for i = 1:6
        diff = abs(params1(i) - params2(i))/mean([params1(i), params2(i)])*100;
        fprintf('  %s: %.4f vs %.4f (差异: %.1f%%)\n', ...
            param_names{i}, params1(i), params2(i), diff);
    end
    
    % 3. 物理合理性检验
    fprintf('\n物理合理性检验:\n');
    
    % 检查厚度范围
    if thickness1 > 1 && thickness1 < 200 && thickness2 > 1 && thickness2 < 200
        fprintf('? 厚度在合理范围内 (1-200 μm)\n');
    else
        fprintf('? 厚度可能超出合理范围\n');
    end
    
    % 检查折射率范围（碳化硅典型值2.5-2.8）
    lambda_test = 2.0;  % 测试波长2μm
    n_test1 = calculate_sellmeier_n(lambda_test, params1);
    n_test2 = calculate_sellmeier_n(lambda_test, params2);
    
    if n_test1 > 2.3 && n_test1 < 3.0 && n_test2 > 2.3 && n_test2 < 3.0
        fprintf('? 折射率在合理范围内 (2.3-3.0)\n');
        fprintf('  λ=2μm时: n1=%.3f, n2=%.3f\n', n_test1, n_test2);
    else
        fprintf('? 折射率可能异常\n');
    end
    
    % 4. 建议
    fprintf('\n建议:\n');
    if thickness_diff_percent < 5
        fprintf('? 结果可靠，两个角度测量结果一致\n');
        fprintf('? 推荐厚度值: %.2f ± %.2f μm\n', ...
            mean([thickness1, thickness2]), ...
            abs(thickness1-thickness2)/2);
    else
        fprintf('? 建议检查实验条件或重新测量\n');
        fprintf('? 可能需要考虑多光束干涉效应\n');
    end
end