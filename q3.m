clc;
clear;
close all;

% 文件配置
file1.path = '附件3.xlsx';
file1.angle = 10;
file2.path = '附件4.xlsx';
file2.angle = 15;

% 全局配置
config.waveN_fit_min = 400; % 拟合的起始波数 (cm^-1)

%% 模型参数与初值设定
config.n1_init = 3.48; % 外延层折射率估算初值

% 外延层: 柯西模型 n(lambda) = A + B / lambda^2
config.cauchyParam_init = [3.42, 0.05];

% 外延层: 吸收模型 k(lambda) = A * lambda^B
config.k1Param_A_init = 1e-5;
config.k1Param_B_init = 2.0;

% 衬底: Drude模型参数 [nu_p (cm^-1), Gamma (cm^-1)]
% nu_p: 等离子体频率，与掺杂浓度有关
% Gamma: 阻尼系数 (散射率)
config.drudeParam_init = [2000, 200];
config.epsilon_inf = 11.7; % 硅的高频介电常数 (固定值)

%% 初值估算
disp('计算厚度初值');
% FFT建议在干涉条纹明显的短波/高波数区域进行
data1 = readmatrix(file1.path);
thk_init1 = fft_thk_estimate(data1(:,1), data1(:,2)/100, config.n1_init, 2000, file1.angle);
data2 = readmatrix(file2.path);
thk_init2 = fft_thk_estimate(data2(:,1), data2(:,2)/100, config.n1_init, 2000, file2.angle);
config.thk_init = mean([thk_init1, thk_init2]);
fprintf('文件3 (10°) FFT估算厚度: %.2f μm\n', thk_init1);
fprintf('文件4 (15°) FFT估算厚度: %.2f μm\n', thk_init2);
fprintf('平均厚度初值: %.2f μm\n', config.thk_init);
fprintf('\n');

%% 数据处理与拟合
disp('--- 开始处理附件3 (10°) ---')
[result1] = process(data1, file1.angle, config);
fprintf('\n');
disp('--- 开始处理附件4 (15°) ---')
[result2] = process(data2, file2.angle, config);

%% 结果分析
analyze_results(result1, result2, file1.angle, file2.angle);

%% 核心处理函数
function [result] = process(data, angle, config)
    waveNum = data(:, 1);
    R = data(:, 2) / 100;
    waveLen_full = 10000 ./ waveNum;

    filter = waveNum > config.waveN_fit_min;
    waveNum_fit = waveNum(filter);
    R_fit = R(filter);
    waveLen_fit = waveLen_full(filter);

    %% 参数向量 [thk, cauchy_A, cauchy_B, k1_A, k1_B, drude_nu_p, drude_Gamma]
    x0 = [config.thk_init, config.cauchyParam_init, config.k1Param_A_init, config.k1Param_B_init, config.drudeParam_init];
    lb = [config.thk_init*0.9, 3.4, 0,    0,    0,   100,  10];
    ub = [config.thk_init*1.1, 3.5, 0.2,  1e-3, 4,   4000, 1000];

    [x_optimal, R_squared_fit] = global_fit_all_parameters(x0, lb, ub, waveNum_fit, R_fit, waveLen_fit, angle, config.epsilon_inf);

    %% 从优化结果中分解参数
    result.thk = x_optimal(1);
    result.cauchyParam = x_optimal(2:3);
    result.k1Param = x_optimal(4:5);
    result.drudeParam = x_optimal(6:7);

    % 计算全谱段的光学常数
    result.n1_complex_full = calculate_n1_complex(waveLen_full, result.cauchyParam, result.k1Param);
    result.n2_complex_full = calculate_n2_complex_drude(waveNum, result.drudeParam, config.epsilon_inf);

    fprintf('   最终厚度: %.2f μm\n', result.thk);
    fprintf('   最终外延层 n_complex (在 6μm): %.3f + %.4fi\n', real(interp1(waveLen_full, result.n1_complex_full, 6)), imag(interp1(waveLen_full, result.n1_complex_full, 6)));
    fprintf('   最终衬底Drude参数 ν_p: %.1f cm??, Γ: %.1f cm??\n', result.drudeParam(1), result.drudeParam(2));
    fprintf('   在拟合区域的拟合优度 R?: %.4f\n', R_squared_fit);

    plot_fit_res(waveNum, R, result.n1_complex_full, result.n2_complex_full, result.thk, angle, config);
    plot_complex_refractive_index(waveLen_full, waveNum, result.n1_complex_full, '外延层');
    plot_complex_refractive_index(waveLen_full, waveNum, result.n2_complex_full, '衬底 (Drude模型)');
end

%% 全局优化函数
function [x_optimal, R_squared] = global_fit_all_parameters(x0, lb, ub, waveNum_fit, R_fit, waveLen_fit, angle, epsilon_inf)
    theta0_rad = angle * pi / 180;
    model_func = @(x, k) model_R(x, k, waveLen_fit, theta0_rad, epsilon_inf);
    options = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxIterations', 1000, ...
        'FunctionTolerance', 1e-9, 'StepTolerance', 1e-10, ...
        'Algorithm', 'trust-region-reflective', 'UseParallel', true);
    [x_optimal, resnorm] = lsqcurvefit(model_func, x0, waveNum_fit, R_fit, lb, ub, options);
    
    SS_tot = sum((R_fit - mean(R_fit)).^2);
    R_squared = 1 - resnorm / SS_tot;
end

%% FFT 厚度估算 (通用函数)
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
    r_win = r_uniform(:) .* window(:);
    fft_result = fft(r_win);
    fft_power = abs(fft_result(1:N/2)).^2;
    dk = mean(diff(k_uniform));
    thk_axis = (0:N/2-1) * 10000 / (2 * n_avg * cos_theta1 * N * dk);
    search_range = find(thk_axis > 5 & thk_axis < 200);
    if isempty(search_range)
        thk = 20; return;
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
    if DatePower2 > 0, f_corrected = DatePower1 / DatePower2;
    else, f_corrected = max_idx_rough; end
    thk = (f_corrected - 1) * 10000 / (2 * n_avg * cos_theta1 * N * dk);
end

%% 结果分析函数
function analyze_results(res1, res2, angle1, angle2)
    disp('--- 最终结果对比 ---');
    fprintf('文件3 (%d°) -> 厚度: %.2f μm, ν_p: %.1f, Γ: %.1f\n', angle1, res1.thk, res1.drudeParam(1), res1.drudeParam(2));
    fprintf('文件4 (%d°) -> 厚度: %.2f μm, ν_p: %.1f, Γ: %.1f\n', angle2, res2.thk, res2.drudeParam(1), res2.drudeParam(2));
    thk1 = res1.thk; thk2 = res2.thk;
    thk_diff_percent = abs(thk1-thk2)/mean([thk1,thk2])*100;
    fprintf('两个角度测得的厚度分别为 %.2f μm 和 %.2f μm，相对差异为 %.2f%%。\n', thk1, thk2, thk_diff_percent);
    if thk_diff_percent < 5
        fprintf('   -> 结论: 厚度一致性良好，结果可靠。\n');
        fprintf('   -> 推荐厚度值: %.2f ± %.2f μm\n', mean([thk1, thk2]), std([thk1, thk2]));
    else
        fprintf('   -> 结论: 厚度差异较大。\n');
    end
end

%% 物理模型与绘图函数

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
    % 外延层: 柯西模型 + 吸收项
    A = cauchyParam(1); B = cauchyParam(2);
    n_real = A + B ./ (waveLen.^2);
    k1_A = k1Param(1); k1_B = k1Param(2);
    k_imag = k1_A * waveLen.^k1_B;
    n1_complex = n_real + 1i * k_imag;
end

function n2_complex = calculate_n2_complex_drude(waveNum, drudeParam, epsilon_inf)
    % 衬底: Drude 模型
    nu_p = drudeParam(1);  % 等离子体频率 (cm^-1)
    Gamma = drudeParam(2); % 阻尼系数 (cm^-1)
    nu = waveNum(:);
    epsilon_complex = epsilon_inf - (nu_p^2) ./ (nu.^2 + 1i * Gamma * nu);
    n2_complex = sqrt(epsilon_complex);
end

function R = model_R_vectorized(thk, n1_complex, n2_complex, waveNum, theta0)
    % 核心物理模型 (通用)
    n0 = 1.0;
    n1 = n1_complex(:);
    n2 = n2_complex(:); % n2 可以是复数
    waveNum = waveNum(:);
    sin_theta1 = n0 * sin(theta0) ./ n1;
    cos_theta1 = sqrt(1 - sin_theta1.^2);
    sin_theta2 = n0 * sin(theta0) ./ n2;
    cos_theta2 = sqrt(1 - sin_theta2.^2); % n2是复数, cos_theta2也将是复数

    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    r12_s = (n1.*cos_theta1 - n2.*cos_theta2) ./ (n1.*cos_theta1 + n2.*cos_theta2);
    r01_p = (n1*cos(theta0) - n0*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    r12_p = (n2.*cos_theta1 - n1.*cos_theta2) ./ (n2.*cos_theta1 + n1.*cos_theta2);

    delta = 4 * pi * n1 .* thk .* cos_theta1 .* waveNum / 10000;
    exp_term = exp(1i*delta);

    % --- 选择干涉模型 ---
    R = calculate_R_double_beam(r01_s, r12_s, r01_p, r12_p, exp_term); % 使用双光束干涉
    % R = calculate_R_multi_beam(r01_s, r12_s, r01_p, r12_p, exp_term);  % 使用多光束干涉
end

%% 干涉模型 (通用函数)
function R = calculate_R_double_beam(r01_s, r12_s, r01_p, r12_p, exp_term)
    % 双光束干涉模型
    r_s_total = r01_s + r12_s .* exp_term;
    r_p_total = r01_p + r12_p .* exp_term;
    R_s = abs(r_s_total).^2;
    R_p = abs(r_p_total).^2;
    R = real((R_s + R_p) / 2);
end

function R = calculate_R_multi_beam(r01_s, r12_s, r01_p, r12_p, exp_term)
    % 多光束干涉模型
    r_s_total = (r01_s + r12_s .* exp_term) ./ (1 + r01_s .* r12_s .* exp_term);
    r_p_total = (r01_p + r12_p .* exp_term) ./ (1 + r01_p .* r12_p .* exp_term);
    R_s = abs(r_s_total).^2;
    R_p = abs(r_p_total).^2;
    R = real((R_s + R_p) / 2);
end

%% 绘图函数 (通用风格)
function plot_complex_refractive_index(waveLen, waveNum, n_complex_full, layer_name)
    figure('Name', ['拟合得到的' layer_name '复折射率色散曲线'], 'Position', [100, 100, 1000, 450]);
    n_real = real(n_complex_full);
    k_imag = imag(n_complex_full);
    t = tiledlayout(1, 2);
    
    % n vs Wavelength/Wavenumber
    ax1 = nexttile;
    colororder(ax1, {'blue', 'red'});
    yyaxis(ax1, 'left');
    plot(waveNum, n_real, 'b-', 'LineWidth', 1.5);
    ylabel(['折射率实部 n (' layer_name ')']);
    ylim([min(n_real)*0.95, max(n_real)*1.05]);
    
    yyaxis(ax1, 'right');
    plot(waveNum, k_imag, 'r-', 'LineWidth', 1.5);
    ylabel(['消光系数 k (' layer_name ')']);
    
    xlabel('波数 (cm^{-1})'); grid on; xlim([min(waveNum), max(waveNum)]);
    title('光学常数 vs 波数');

    % n vs Wavelength
    ax2 = nexttile;
    colororder(ax2, {'blue', 'red'});
    yyaxis(ax2, 'left');
    plot(waveLen, n_real, 'b-', 'LineWidth', 1.5);
    ylabel(['折射率实部 n (' layer_name ')']);
    ylim([min(n_real)*0.95, max(n_real)*1.05]);

    yyaxis(ax2, 'right');
    plot(waveLen, k_imag, 'r-', 'LineWidth', 1.5);
    ylabel(['消光系数 k (' layer_name ')']);
    
    xlabel('波长 (μm)'); grid on;
    title('光学常数 vs 波长');
end

function plot_fit_res(waveNum, R, n1_complex_full, n2_complex_full, thk, angle, config)
    theta0 = angle * pi / 180;
    R_fit_full = model_R_vectorized(thk, n1_complex_full, n2_complex_full, waveNum, theta0);
    R_squared_full = 1 - sum((R - R_fit_full).^2) / sum((R - mean(R)).^2);
    
    figure('Name', ['拟合结果分析 (入射角 ' num2str(angle) '°)'], 'Position', [150, 150, 1200, 600]);
    plot(waveNum, R*100, '.', 'Color', [0 .4 .8], 'MarkerSize', 6, 'DisplayName', '实测数据');
    hold on;
    plot(waveNum, R_fit_full*100, 'r-', 'LineWidth', 1.5, 'DisplayName', '拟合模型');
    
    ylim_vals = ylim;
    h = fill([config.waveN_fit_min, max(waveNum), max(waveNum), config.waveN_fit_min], [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
             'k', 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'DisplayName', '拟合区域');
    uistack(h, 'bottom');

    xlabel('波数 (cm^{-1})'); ylabel('反射率 (%)');
    title_str = sprintf('反射率拟合 (θ=%d°), 全局 R?=%.4f', angle, R_squared_full);
    title(title_str);
    legend('Location', 'best'); grid on; xlim([min(waveNum), max(waveNum)]);
end