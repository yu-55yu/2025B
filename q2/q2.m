clc;
clear;
close all;

% --- 文件和配置 ---
file1.path = '附件1.xlsx';
file1.angle = 10;
file2.path = '附件2.xlsx';
file2.angle = 15;

config.waveN_fit_min = 1200;
config.n1_init = 2.58;
config.n2_real_init = 2.55; % 衬底折射率 n2 拟合初值 (现在是实数)
config.k1_A_init = 0.001;   % 外延层吸收模型参数 k1 = A * lambda^B
config.k1_B_init = 2.0;
config.selParam_init = [5.5, 0.2, 0.05, 0.027, 100, 0.01];

% --- 1. 厚度初值估算 ---
disp('计算厚度初值');
data1 = readmatrix(file1.path);
thk_init1 = fft_thk_estimate(data1(:,1), data1(:,2)/100, config.n1_init, config.waveN_fit_min, file1.angle);
data2 = readmatrix(file2.path);
thk_init2 = fft_thk_estimate(data2(:,1), data2(:,2)/100, config.n1_init, config.waveN_fit_min, file2.angle);
config.thk_init = mean([thk_init1, thk_init2]);
fprintf('文件1 (10°) FFT估算厚度: %.2f μm\n', thk_init1);
fprintf('文件2 (15°) FFT估算厚度: %.2f μm\n', thk_init2);
fprintf('平均厚度初值: %.2f μm\n', config.thk_init);

% --- 2. 对每个文件进行拟合处理 ---
fprintf('\n');
disp('--- 开始处理附件1 (10°)')
[result1] = process_file(data1, file1.angle, config);
fprintf('\n');
disp('--- 开始处理附件2 (15°)')
[result2] = process_file(data2, file2.angle, config);

% --- 3. 最终结果对比分析 ---
analyze_results(result1, result2, file1.angle, file2.angle);

%% 核心处理函数
function [result] = process_file(data, incident_angle, config)
    waveNum = data(:, 1);
    R = data(:, 2) / 100;
    waveLen_full = 10000 ./ waveNum;
    filter = waveNum > config.waveN_fit_min;
    waveNum_fit = waveNum(filter);
    R_fit = R(filter);
    waveLen_fit = waveLen_full(filter);

    % thk, n2_r, selParam(6), k1_A, k1_B
    x0 = [config.thk_init, config.n2_real_init, config.selParam_init, config.k1_A_init, config.k1_B_init];
    
    lb = [config.thk_init*0.8, 2.0, 0.1, 0.001, 0.0001, 0.0001, 0.1,  0.001,  0, 0];
    ub = [config.thk_init*1.2, 3.5, 20,  10,     5,       2,       150,  20,     0.01, 4];

    [x_optimal, R_squared_fit] = global_fit_all_parameters(x0, lb, ub, waveNum_fit, R_fit, waveLen_fit, incident_angle);

    result.thk = x_optimal(1);
    result.n2 = x_optimal(2);
    result.selParam = x_optimal(3:8);
    result.k1Param = x_optimal(9:10);
    
    result.n1_complex_full = calculate_n1_complex(waveLen_full, result.selParam, result.k1Param);

    fprintf('  最终厚度: %.2f μm\n', result.thk);
    fprintf('  最终外延层 ñ1 (在 6μm): %.3f + %.4fi\n', real(interp1(waveLen_full, result.n1_complex_full, 6)), imag(interp1(waveLen_full, result.n1_complex_full, 6)));
    fprintf('  最终衬底 n2: %.3f\n', result.n2);
    fprintf('  在拟合区域的拟合优度 R²: %.4f\n', R_squared_fit);

    plot_fitting_results(waveNum, R, result.n1_complex_full, result.thk, result.n2, incident_angle, config.waveN_fit_min);
    plot_refractive_index(waveLen_full, waveNum, result.n1_complex_full);
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
    if nargin < 5
        theta0_deg = 0; 
    end 
        theta0_rad = theta0_deg * pi / 180; 
        cos_theta1 = real(sqrt(1 - (sin(theta0_rad) / n_avg)^2)); 
        filter = waveNum > waveN_fit_min; waveNum_fft = waveNum(filter); 
        R_fft = R(filter); 
        R_ac = R_fft - mean(R_fft); 
        N = 2^nextpow2(8*length(waveNum_fft)); 
        k_uniform = linspace(min(waveNum_fft), max(waveNum_fft), N); 
        r_uniform = interp1(waveNum_fft, R_ac, k_uniform, 'pchip', 'extrap'); 
        window = hann(N); 
        r_windowed = r_uniform(:) .* window(:); fft_result = fft(r_windowed); 
        fft_power = abs(fft_result(1:N/2)).^2; dk = mean(diff(k_uniform)); 
        thk_axis = (0:N/2-1) * 10000 / (2 * n_avg * cos_theta1 * N * dk); 
        search_range = find(thk_axis > 5 & thk_axis < 200); 
    if isempty(search_range)
        thk = 20; 
        return; 
    end
        [~, max_idx_in_range] = max(fft_power(search_range)); 
        max_idx_rough = search_range(max_idx_in_range); 
        correctNum = 3; 
        DatePower1 = 0; 
        DatePower2 = 0; 
        for i = -correctNum:correctNum, idx = max_idx_rough + i; 
            if idx >= 1 && idx <= length(fft_power)
                power = fft_power(idx); 
                DatePower1 = DatePower1 + idx * power; 
                DatePower2 = DatePower2 + power; 
            end
        end
        if DatePower2 > 0
            f_corrected = DatePower1 / DatePower2; 
        else
            f_corrected = max_idx_rough; 
        end
        thk = (f_corrected - 1) * 10000 / (2 * n_avg * cos_theta1 * N * dk);
end

%% 结果分析函数
function analyze_results(res1, res2, angle1, angle2)
    disp('--- 最终结果对比 ---');
    fprintf('文件1 (%d°) -> 厚度: %.2f μm, 衬底 n2: %.3f\n', angle1, res1.thk, res1.n2);
    fprintf('文件2 (%d°) -> 厚度: %.2f μm, 衬底 n2: %.3f\n', angle2, res2.thk, res2.n2);
    thk1 = res1.thk; thk2 = res2.thk;
    thk_diff_percent = abs(thk1-thk2)/mean([thk1,thk2])*100;
    fprintf('两个角度测得的厚度分别为 %.2f μm 和 %.2f μm，相对差异为 %.2f%%。\n', thk1, thk2, thk_diff_percent);
    if thk_diff_percent < 5, fprintf('  -> 结论: 厚度一致性良好，结果可靠。\n'); fprintf('  -> 推荐厚度值: %.2f ± %.2f μm\n', mean([thk1, thk2]), std([thk1, thk2])); else, fprintf('  -> 结论: 厚度差异较大。\n'); end
end

%% --- 物理模型与绘图函数 ---

function R = model_R(params, ~, waveLen, theta0) 
    thk = params(1);
    n2 = params(2);
    selParam = params(3:8);
    k1Param = params(9:10);
    
    n1_complex = calculate_n1_complex(waveLen, selParam, k1Param);
    waveNum_model = 10000 ./ waveLen;
    
    R = model_R_vectorized(thk, n1_complex, n2, waveNum_model, theta0);
end

function n1_complex = calculate_n1_complex(waveLen, selParam, k1Param)
    B1=selParam(1); B2=selParam(2); B3=selParam(3);
    C1=selParam(4); C2=selParam(5); C3=selParam(6);
    lambda_sq = waveLen.^2; eps = 1e-10;
    term1 = B1*lambda_sq./(lambda_sq-C1+eps);
    term2 = B2*lambda_sq./(lambda_sq-C2+eps);
    term3 = B3*lambda_sq./(lambda_sq-C3+eps);
    n_squared = 1 + term1 + term2 + term3;
    n_squared(n_squared < 1) = 1;
    n_real = real(sqrt(n_squared)); 
    n_real(n_real < 1) = 1;
    
    k1_A = k1Param(1);
    k1_B = k1Param(2);
    k_imag = k1_A * waveLen.^k1_B;
    
    n1_complex = n_real + 1i * k_imag;
end

function R = model_R_vectorized(thk, n1_complex, n2_real, waveNum, theta0)
    n0 = 1.0;
    n1 = n1_complex(:);
    n2 = n2_real(:);
    waveNum = waveNum(:);

    sin_theta1 = n0 * sin(theta0) ./ n1;
    cos_theta1 = sqrt(1 - sin_theta1.^2);
    sin_theta2 = n0 * sin(theta0) ./ n2;
    cos_theta2 = sqrt(1 - sin_theta2.^2);

    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    r12_s = (n1.*cos_theta1 - n2.*cos_theta2) ./ (n1.*cos_theta1 + n2.*cos_theta2);
    r01_p = (n1*cos(theta0) - n0*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    r12_p = (n2.*cos_theta1 - n1.*cos_theta2) ./ (n2.*cos_theta1 + n1.*cos_theta2);

    delta = 4 * pi * n1 .* thk .* cos_theta1 .* waveNum / 10000;
    exp_term = exp(1i*delta);

    % 双光束干涉
    r_s_total = r01_s + r12_s .* exp_term;
    r_p_total = r01_p + r12_p .* exp_term;

    R_s = abs(r_s_total).^2;
    R_p = abs(r_p_total).^2;
    R = (R_s + R_p) / 2;
    R = real(R(:));
end


% 绘图函数: 绘制外延层 ñ1 的n和k曲线
function plot_refractive_index(waveLen, waveNum, n1_complex_all)
    % --- 定义颜色 ---
    color_blue = [0, 176, 232] / 255; % #00b0e8
    color_red = [254, 100, 108] / 255; % #fe646c
    
    figure('Name', '拟合得到的外延层复折射率色散曲线', 'Position', [100, 100, 1200, 500]);
    n_real = real(n1_complex_all);
    k_imag = imag(n1_complex_all);
    
    disp('  分析外延层色散特性...');
    dn_dlambda = diff(n_real) ./ diff(waveLen);
    waveLen_mid = (waveLen(1:end-1) + waveLen(2:end)) / 2;
    is_anomalous = dn_dlambda > 0;
    
    if any(is_anomalous)
        start_indices = find(diff([0; is_anomalous]) == 1);
        end_indices = find(diff([is_anomalous; 0]) == -1);
        fprintf('  识别到 %d 个反常色散区域 (Anomalous Dispersion, dn/dλ > 0):\n', length(start_indices));
        for i = 1:length(start_indices)
            fprintf('    - 范围: %.2f μm 到 %.2f μm\n', waveLen_mid(start_indices(i)), waveLen_mid(end_indices(i)));
        end
    else
        disp('  在分析范围内未发现明显的反常色散区。');
    end
    
    % --- n vs 波长 ---
    subplot(2,2,1); 
    plot(waveLen, n_real, '-', 'Color', color_blue, 'LineWidth', 1.5, 'DisplayName', '正常色散区 (dn/dλ < 0)'); 
    hold on;
    n_mid = (n_real(1:end-1) + n_real(2:end)) / 2;
    if any(is_anomalous)
        plot(waveLen_mid(is_anomalous), n_mid(is_anomalous), '.', 'Color', color_red, 'MarkerSize', 8, 'DisplayName', '反常色散区 (dn/dλ > 0)');
    end
    xlabel('波长 (μm)'); ylabel('折射率实部 n_1'); 
    grid on; legend('Location', 'best');
    
    % --- k vs 波长 ---
    subplot(2,2,2); 
    plot(waveLen, k_imag, '-', 'Color', color_red, 'LineWidth', 1.5); 
    xlabel('波长 (μm)'); ylabel('消光系数 k_1'); 
    grid on; set(gca, 'YScale', 'log');
    
    % --- n vs 波数 ---
    subplot(2,2,3); 
    [waveNum_sorted, sort_idx] = sort(waveNum);
    n_real_sorted = n_real(sort_idx);
    dn_dk = diff(n_real_sorted) ./ diff(waveNum_sorted);
    waveNum_mid = (waveNum_sorted(1:end-1) + waveNum_sorted(2:end)) / 2;
    n_mid_k = (n_real_sorted(1:end-1) + n_real_sorted(2:end)) / 2;
    is_anomalous_k = dn_dk < 0;
    
    plot(waveNum, n_real, '-', 'Color', color_blue, 'LineWidth', 1.5, 'DisplayName', '正常色散区 (dn/dk > 0)');
    hold on;
    if any(is_anomalous_k)
        plot(waveNum_mid(is_anomalous_k), n_mid_k(is_anomalous_k), '.', 'Color', color_red, 'MarkerSize', 8, 'DisplayName', '反常色散区 (dn/dk < 0)');
    end
    xlabel('波数 (cm^{-1})'); ylabel('折射率实部 n_1'); 
    grid on; xlim([min(waveNum), max(waveNum)]); legend('Location', 'best');

    % --- k vs 波数 ---
    subplot(2,2,4); 
    plot(waveNum, k_imag, '-', 'Color', color_red, 'LineWidth', 1.5); 
    xlabel('波数 (cm^{-1})'); ylabel('消光系数 k_1'); 
    grid on; xlim([min(waveNum), max(waveNum)]); set(gca, 'YScale', 'log');
end

% 绘图函数: 绘制拟合结果 (*** 已修改 ***)
function plot_fitting_results(waveNum, R, n1_complex_full, thk, n2_real, incident_angle, waveN_fit_min)
    % --- 定义颜色 ---
    color_blue = [0, 176, 232] / 255; % #00b0e8
    color_red = [254, 100, 108] / 255; % #fe646c

    theta0 = incident_angle * pi / 180;
    R_fitted = model_R_vectorized(thk, n1_complex_full, n2_real, waveNum, theta0);
    
    % --- 图 1: 全局拟合结果 ---
    figure('Name', ['拟合结果分析 (入射角 ' num2str(incident_angle) '°)'], 'Position', [150, 150, 1200, 600]);
    plot(waveNum, R*100, '.', 'Color', color_blue, 'MarkerSize', 5, 'DisplayName', '实测数据');
    hold on;
    plot(waveNum, R_fitted*100, '-', 'Color', color_red, 'LineWidth', 1.5, 'DisplayName', '拟合模型');
    ylim_vals = ylim;
    h = fill([waveN_fit_min, max(waveNum), max(waveNum), waveN_fit_min], [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], 'k', 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'DisplayName', '拟合区域');
    uistack(h, 'bottom');
    xlabel('波数 (cm^{-1})'); ylabel('反射率 (%)');
    title(['拟合结果 (入射角 ' num2str(incident_angle) '°)']);
    legend('Location', 'best'); grid on; xlim([min(waveNum), max(waveNum)]);

    % --- 新增 图 2: 1500-2000 cm-1 局部放大图 ---
    figure('Name', ['局部放大 (入射角 ' num2str(incident_angle) '°)'], 'Position', [200, 200, 900, 500]);
    plot(waveNum, R*100, '.', 'Color', color_blue, 'MarkerSize', 8, 'DisplayName', '实测数据');
    hold on;
    plot(waveNum, R_fitted*100, '-', 'Color', color_red, 'LineWidth', 1.5, 'DisplayName', '拟合模型');
    xlabel('波数 (cm^{-1})'); ylabel('反射率 (%)');
    title(['局部放大图 (1500-2000 cm^{-1})']);
    legend('Location', 'best'); grid on;
    xlim([1500, 2000]); % 设置X轴范围
end