clc;
clear;
close all;

% --- 1. 修改文件路径以适用于附件3和附件4 (硅) ---
file1.path = '附件3.xlsx';
file1.angle = 10;
file2.path = '附件4.xlsx';
file2.angle = 15;

% --- 2. 根据您提供的图片信息，更新配置参数 ---
config.waveN_fit_min = 1000;      % 硅的数据在较低波数仍然清晰，可适当放宽
% [cite: 1]
config.n1_init = 3.43;            % 用于FFT估算的硅外延层平均折射率 (取1.5-5μm范围的推荐值)
% [cite: 1]
config.n2_real_init = 3.43;         % 衬底折射率实部 n2 拟合初值
config.k2_imag_init = 10;         % 硅衬底吸收相对较弱，给一个较小的初值
% [cite: 1]
% 根据图片提供的硅Sellmeier公式 n^2-1 = B1*L^2/(L^2-C1_sqrt^2) + ... 设置初值
% 其中 B对应分子系数, C对应分母波长的平方 (C = lambda^2)
B1 = 10.6684; B2 = 0.00304; B3 = 1.5413; % [cite: 1]
C1_sqrt = 0.3015; C2_sqrt = 1.1348; C3_sqrt = 1104; % [cite: 1]
config.selParam_init = [B1, B2, B3, C1_sqrt^2, C2_sqrt^2, C3_sqrt^2];


disp('计算厚度初值');
data1 = readmatrix(file1.path);
thk_init1 = fft_thk_estimate(data1(:,1), data1(:,2)/100, config.n1_init, config.waveN_fit_min, file1.angle);

data2 = readmatrix(file2.path);
thk_init2 = fft_thk_estimate(data2(:,1), data2(:,2)/100, config.n1_init, config.waveN_fit_min, file2.angle);

config.thk_init = mean([thk_init1, thk_init2]);
fprintf('文件3 (10°) FFT估算厚度: %.2f μm\n', thk_init1);
fprintf('文件4 (15°) FFT估算厚度: %.2f μm\n', thk_init2);
fprintf('平均厚度初值: %.2f μm\n', config.thk_init);

fprintf('\n');
disp('--- 开始处理附件3 (10°) ---')
[result1.thk, result1.n2_complex, result1.params] = process_file(data1, file1.angle, config);

fprintf('\n');
disp('--- 开始处理附件4 (15°) ---')
[result2.thk, result2.n2_complex, result2.params] = process_file(data2, file2.angle, config);

analyze_results(result1, result2, file1.angle, file2.angle);

%% 核心处理函数
function [thk, n2_complex, selParam] = process_file(data, incident_angle, config)
    waveNum = data(:, 1);  % 全波段波数
    R = data(:, 2) / 100;  % 全波段反射率
    waveLen_full = 10000 ./ waveNum;

    filter = waveNum > config.waveN_fit_min;
    waveNum_fit = waveNum(filter);
    R_fit = R(filter);
    waveLen_fit = waveLen_full(filter);

    % 参数顺序: [厚度, n2实部, k2虚部, 6个Sellmeier系数]
    x0 = [config.thk_init, config.n2_real_init, config.k2_imag_init, config.selParam_init];
    
    % 参数边界 (为硅材料调整)
    lb = [config.thk_init*0.8, 3.0, 10,    10, 0.001, 1,   0.09, 1.2, 1e5]; % Lower bounds
    ub = [config.thk_init*1.2, 4.0, 30,  11, 0.005, 2,   0.1,  1.3, 2e6]; % Upper bounds

    [x_optimal, R_squared_fit] = global_fit_all_parameters(x0, lb, ub, waveNum_fit, R_fit, waveLen_fit, incident_angle);

    thk = x_optimal(1);
    n2_complex = x_optimal(2) + 1i * x_optimal(3); % 组合成复数折射率
    selParam = x_optimal(4:9);

    fprintf('  最终厚度: %.2f μm\n', thk);
    fprintf('  最终衬底复折射率 ?2: %.3f + %.4fi\n', real(n2_complex), imag(n2_complex));
    fprintf('  在拟合区域的拟合优度 R?: %.4f\n', R_squared_fit);

    n1_all = calculate_sellmeier_n(waveLen_full, selParam);
    plot_fitting_results(waveNum, R, n1_all, thk, n2_complex, incident_angle, config.waveN_fit_min);
    plot_refractive_index(waveLen_full, waveNum, n1_all);
end

%% 全局优化函数
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

%% FFT 厚度估算
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

%% 结果分析函数
function analyze_results(res1, res2, angle1, angle2)
    disp('--- 最终结果对比 ---');
    fprintf('文件3 (%d°) -> 厚度: %.2f μm, 衬底 ?2: %.3f + %.4fi\n', angle1, res1.thk, real(res1.n2_complex), imag(res1.n2_complex));
    fprintf('文件4 (%d°) -> 厚度: %.2f μm, 衬底 ?2: %.3f + %.4fi\n', angle2, res2.thk, real(res2.n2_complex), imag(res2.n2_complex));

    thk1 = res1.thk; thk2 = res2.thk;
    thk_diff_percent = abs(thk1-thk2)/mean([thk1,thk2])*100;
    fprintf('两个角度测得的厚度分别为 %.2f μm 和 %.2f μm，相对差异为 %.2f%%。\n', thk1, thk2, thk_diff_percent);
    
    if thk_diff_percent < 2
        fprintf('  -> 结论: 厚度一致性良好，结果可靠。\n');
        fprintf('  -> 推荐厚度值: %.2f ± %.2f μm\n', mean([thk1, thk2]), std([thk1, thk2]));
    else
        fprintf('  -> 结论: 厚度差异较大，建议检查模型或数据。\n');
    end
end

%% --- 物理模型与绘图函数 ---

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
    figure('Name', '拟合得到的外延层折射率色散曲线', 'Position', [100, 100, 1200, 500]);
    subplot(1,2,1); plot(waveLen, n1_all, 'b-', 'LineWidth', 2); xlabel('波长 (μm)'); ylabel('外延层折射率 n_1'); title('折射率 vs 波长'); grid on;
    subplot(1,2,2); plot(waveNum, n1_all, 'r-', 'LineWidth', 2); xlabel('波数 (cm^{-1})'); ylabel('外延层折射率 n_1'); title('折射率 vs 波数'); grid on; xlim([min(waveNum), max(waveNum)]);
end

function plot_fitting_results(waveNum, R, n1, thk, n2_complex, incident_angle, waveN_fit_min)
    theta0 = incident_angle * pi / 180;
    R_fitted = model_R_vectorized(thk, n1, n2_complex, waveNum, theta0);
    R_squared_full = 1 - sum((R - R_fitted).^2) / sum((R - mean(R)).^2);

    figure('Name', ['拟合结果分析 (入射角 ' num2str(incident_angle) '°)'], 'Position', [150, 150, 1200, 600]);
    plot(waveNum, R*100, 'b.', 'MarkerSize', 5, 'DisplayName', '实测数据');
    hold on;
    plot(waveNum, R_fitted*100, 'r-', 'LineWidth', 2, 'DisplayName', '拟合模型');
    ylim_vals = ylim;
    h = fill([waveN_fit_min, max(waveNum), max(waveNum), waveN_fit_min], [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], 'k', 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'DisplayName', '拟合区域');
    uistack(h, 'bottom');
    xlabel('波数 (cm^{-1})'); ylabel('反射率 (%)');
    title_str = sprintf('反射率拟合结果 (θ=%d°, d=%.2f μm, ?_2=%.3f+%.4fi, 全局R?=%.4f)', ...
        incident_angle, thk, real(n2_complex), imag(n2_complex), R_squared_full);
    title(title_str);
    legend('Location', 'best'); grid on; xlim([min(waveNum), max(waveNum)]);
end