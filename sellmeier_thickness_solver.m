
function sellmeier_thickness_solver()
    clc;
    clear;
    close all;
    
    % 处理附件    
    [thickness1, n2_1, params1, ~] = process_file('附件1.xlsx', 10);    
    [thickness2, n2_2, params2, ~] = process_file('附件2.xlsx', 15);
    
    % 结果
    fprintf('\n==================================================\n');
    fprintf('最终结果对比分析\n');
    fprintf('==================================================\n');
    fprintf('附件1 -> 厚度: %.2f μm, 衬底折射率 n2: %.3f\n', thickness1, n2_1);
    fprintf('附件2 -> 厚度: %.2f μm, 衬底折射率 n2: %.3f\n', thickness2, n2_2);
    fprintf('厚度差异: %.2f μm (%.2f%%)\n', ...
        abs(thickness1-thickness2), ...
        abs(thickness1-thickness2)/mean([thickness1,thickness2])*100);
    
    % 结果可靠性分析
    analyze_reliability(thickness1, thickness2, params1, params2);
end


function [thickness, n2, sellmeier_params, n1_all] = process_file(filename, incident_angle)
    data = readmatrix(filename);
    wavenumber_full = data(:, 1);  % 全波段波数
    reflectance_full = data(:, 2) / 100;  % 全波段反射率
    wavelength_full = 10000 ./ wavenumber_full;  % 全波段波长
    

    %筛选
    fit_mask = wavenumber_full > 1500;
    wavenumber_fit = wavenumber_full(fit_mask);
    reflectance_fit = reflectance_full(fit_mask);
    wavelength_fit = wavelength_full(fit_mask);
    fprintf('  已筛选波数 > 1500 cm^-1 的数据, 共 %d 点用于处理。\n', length(wavenumber_fit));

    fprintf('\n步骤3: 估算参数初始值...\n');
    n_avg_sic = 2.5835; 

    thickness_init = fft_thickness_estimate_improved(wavenumber_fit, reflectance_fit, n_avg_sic);
    fprintf('  FFT 估计厚度初值 (基于高波数数据): d_init = %.2f μm\n', thickness_init);
    
    n2_init = 2.55;
    fprintf('  衬底折射率初值: n2_init = %.2f\n', n2_init);
    
    B1_init = 5.5394;  C1_init = 0.026945;
    B2_init = 0.20;  C2_init = 100;
    B3_init = 0.05;  C3_init = 0.01;
    sellmeier_params_init = [B1_init, B2_init, B3_init, C1_init, C2_init, C3_init];
    fprintf('  塞尔迈耶尔参数初值已设置 (基于文献)\n');

    fprintf('\n步骤4: 执行一体化全局拟合...\n');
    
    x0 = [thickness_init, n2_init, sellmeier_params_init];
    lb = [thickness_init*0.8, 2.0,  0.1, 0.001, 0.0001, 0.0001, 0.1,  0.001];
    ub = [thickness_init*1.2, 3.5,  20,  10,    5,      2,      150,  20];
    

    [x_optimal, R_squared_fit] = global_fit_all_parameters(x0, lb, ub, wavenumber_fit, reflectance_fit, wavelength_fit, incident_angle);
    

    thickness = x_optimal(1);
    n2 = x_optimal(2);
    sellmeier_params = x_optimal(3:8);
    
    fprintf('\n--- 拟合完成 ---\n');
    fprintf('最终厚度: %.2f μm\n', thickness);
    fprintf('最终衬底折射率 n2: %.3f\n', n2);
    fprintf('最终塞尔迈耶尔参数:\n');
    fprintf('  B1=%.4f, B2=%.4f, B3=%.4f\n', sellmeier_params(1), sellmeier_params(2), sellmeier_params(3));
    fprintf('  C1=%.4f, C2=%.4f, C3=%.4f\n', sellmeier_params(4), sellmeier_params(5), sellmeier_params(6));
    fprintf('在拟合区域(>1500 cm^-1)的拟合优度 R²: %.4f\n', R_squared_fit);

    fprintf('\n步骤5: 计算最终折射率并绘制全波段对比图...\n');
    

    n1_all = calculate_sellmeier_n(wavelength_full, sellmeier_params);
    
    plot_refractive_index(wavelength_full, wavenumber_full, n1_all);
    plot_fitting_results_improved(wavenumber_full, reflectance_full, n1_all, thickness, n2, incident_angle);
end

%  全局优化函数
function [x_optimal, R_squared] = global_fit_all_parameters(x0, lb, ub, wavenumber_fit, reflectance_fit, wavelength_fit, incident_angle)
    
    theta0_rad = incident_angle * pi / 180;
    model_func = @(x, k) model_reflectance_full(x, k, wavelength_fit, theta0_rad);
    
    options = optimoptions('lsqcurvefit', ...
        'Display', 'iter', 'MaxIterations', 1000, 'FunctionTolerance', 1e-9, ...
        'StepTolerance', 1e-10, 'Algorithm', 'trust-region-reflective', 'UseParallel', true);
        
    [x_optimal, ~] = lsqcurvefit(model_func, x0, wavenumber_fit, reflectance_fit, lb, ub, options);
    
    R_fitted = model_func(x_optimal, wavenumber_fit);
    SS_tot = sum((reflectance_fit - mean(reflectance_fit)).^2);
    SS_res = sum((reflectance_fit - R_fitted).^2);
    R_squared = 1 - SS_res / SS_tot;
end


%  物理模型与辅助函数
function R = model_reflectance_full(params, wavenumber, wavelength, theta0)
    thickness = params(1);
    n2 = params(2);
    sellmeier_params = params(3:8);
    n1 = calculate_sellmeier_n(wavelength, sellmeier_params);
    R = model_reflectance_vectorized([thickness, n2], wavenumber, n1, theta0);
end

function R = model_reflectance_vectorized(params, wavenumber, n1, theta0)
    thickness = params(1); n2 = params(2); n0 = 1.0;
    n1 = n1(:); wavenumber = wavenumber(:);
    sin_theta1 = n0 * sin(theta0) ./ n1; cos_theta1 = real(sqrt(1 - sin_theta1.^2));
    sin_theta2 = n0 * sin(theta0) / n2; cos_theta2 = real(sqrt(1 - sin_theta2^2));
    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    r12_s = (n1.*cos_theta1 - n2*cos_theta2) ./ (n1.*cos_theta1 + n2*cos_theta2);
    r01_p = (n1*cos(theta0) - n0.*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    r12_p = (n2*cos_theta1 - n1.*cos_theta2) ./ (n2*cos_theta1 + n1.*cos_theta2);
    delta = 4 * pi * n1 .* thickness .* cos_theta1 .* wavenumber / 10000;
    R_s = abs((r01_s + r12_s .* exp(1i*delta)) ./ (1 + r01_s .* r12_s .* exp(1i*delta))).^2;
    R_p = abs((r01_p + r12_p .* exp(1i*delta)) ./ (1 + r01_p .* r12_p .* exp(1i*delta))).^2;
    R = (R_s + R_p) / 2; R = real(R(:));
end

function n = calculate_sellmeier_n(wavelength, params)
    B1 = params(1); B2 = params(2); B3 = params(3);
    C1 = params(4); C2 = params(5); C3 = params(6);
    lambda_sq = wavelength.^2; eps = 1e-10;
    term1 = B1 * lambda_sq ./ (lambda_sq - C1 + eps);
    term2 = B2 * lambda_sq ./ (lambda_sq - C2 + eps);  
    term3 = B3 * lambda_sq ./ (lambda_sq - C3 + eps);
    n_squared = 1 + term1 + term2 + term3;
    n_squared(n_squared < 1) = 1; 
    n = real(sqrt(n_squared)); n(n < 1) = 1;
end

function thickness = fft_thickness_estimate_improved(wavenumber, reflectance, n_avg)
    reflectance_ac = reflectance - mean(reflectance);
    N = 2^nextpow2(8*length(wavenumber));
    k_uniform = linspace(min(wavenumber), max(wavenumber), N);
    r_uniform = interp1(wavenumber, reflectance_ac, k_uniform, 'pchip', 'extrap');
    window = hann(N); r_windowed = r_uniform(:) .* window(:);
    fft_result = fft(r_windowed); fft_power = abs(fft_result(1:N/2)).^2;
    dk = mean(diff(k_uniform));
    thickness_axis = (0:N/2-1) * 10000 / (2 * n_avg * N * dk);
    search_range = find(thickness_axis > 5 & thickness_axis < 200);
    if isempty(search_range), thickness = 20; warning('FFT未能找到明显峰值，使用默认厚度初值20um'); return; end
    [~, max_idx] = max(fft_power(search_range));
    thickness = thickness_axis(search_range(max_idx));
end

function plot_refractive_index(wavelength, wavenumber, n1_all)
    figure('Name', '拟合得到的折射率色散曲线', 'Position', [100, 100, 1200, 500]);
    subplot(1,2,1); plot(wavelength, n1_all, 'b-', 'LineWidth', 2); xlabel('波长 (μm)'); ylabel('外延层折射率 n_1'); title('折射率 vs 波长'); grid on;
    subplot(1,2,2); plot(wavenumber, n1_all, 'r-', 'LineWidth', 2); xlabel('波数 (cm^{-1})'); ylabel('外延层折射率 n_1'); title('折射率 vs 波数'); grid on; xlim([min(wavenumber), max(wavenumber)]);
end

function plot_fitting_results_improved(wavenumber, reflectance, n1, thickness, n2, incident_angle)
    theta0 = incident_angle * pi / 180;
    R_fitted = model_reflectance_vectorized([thickness, n2], wavenumber, n1, theta0);
    R_squared_full = 1 - sum((reflectance - R_fitted).^2) / sum((reflectance - mean(reflectance)).^2);
    
    figure('Name', '全局拟合结果分析 (高波数区域)', 'Position', [150, 150, 1200, 600]);
    plot(wavenumber, reflectance*100, 'b.', 'MarkerSize', 5, 'DisplayName', '实测数据 (全波段)');
    hold on;
    plot(wavenumber, R_fitted*100, 'r-', 'LineWidth', 2, 'DisplayName', '拟合模型 (外推至全波段)');
    
    ylim_vals = ylim;
    h = fill([1500, max(wavenumber), max(wavenumber), 1500], [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], 'k', 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'DisplayName', '拟合与处理区域');
    uistack(h, 'bottom');
    
    xlabel('波数 (cm^{-1})'); ylabel('反射率 (%)');
    title(sprintf('反射率拟合结果 (θ=%d°, d=%.2f μm, n_2=%.3f, 全局R²=%.4f)', ...
        incident_angle, thickness, n2, R_squared_full));
    legend('Location', 'best'); grid on; xlim([min(wavenumber), max(wavenumber)]);
end

function analyze_reliability(thickness1, thickness2, params1, params2)
    fprintf('\n=== 结果可靠性分析 ===\n');
    thickness_diff_percent = abs(thickness1-thickness2)/mean([thickness1,thickness2])*100;
    fprintf('厚度一致性: 两个角度测得的厚度分别为 %.2f μm 和 %.2f μm，相对差异为 %.2f%%。\n', thickness1, thickness2, thickness_diff_percent);
    if thickness_diff_percent < 5, fprintf('  -> 结论: 一致性良好，结果可靠。\n');
        fprintf('  -> 推荐厚度值: %.2f ± %.2f μm\n', mean([thickness1, thickness2]), abs(thickness1-thickness2)/2);
    else, fprintf('  -> 结论: 差异较大，建议检查模型或数据。\n'); end
    fprintf('\n塞尔迈耶尔参数稳定性:\n');
    param_names = {'B1', 'B2', 'B3', 'C1', 'C2', 'C3'}; is_stable = true;
    for i = 1:6, diff = abs(params1(i) - params2(i))/mean([params1(i), params2(i)])*100;
        if diff > 50, is_stable = false; end
        fprintf('  %s: %.4f vs %.4f (差异: %.1f%%)\n', param_names{i}, params1(i), params2(i), diff); end
    if is_stable, fprintf('  -> 结论: 参数在不同角度下拟合结果基本稳定。\n');
    else, fprintf('  -> 结论: 部分参数在不同角度下拟合结果差异较大，但只要最终的n(λ)曲线和厚度d一致即可。\n'); end
end
