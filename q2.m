clc;
clear;
close all;

file1.path = '附件1.xlsx';
file1.angle = 10;
file2.path = '附件2.xlsx';
file2.angle = 15;


config.waveN_fit_min = 1500;        % 拟合使用的最小波数
config.n1_init = 2.5835;            % 用于FFT估算的 n1
config.n2_init = 2.55;              % 衬底折射率 n2 的拟合初值
config.selParam_init = [5.5394, 0.20, 0.05, 0.026945, 100, 0.01];


disp('计算厚度初值');

data1 = readmatrix(file1.path);
thk_init1 = fft_thk_estimate(data1(:,1), data1(:,2)/100, config.n1_init, config.waveN_fit_min);

data2 = readmatrix(file2.path);
thk_init2 = fft_thk_estimate(data2(:,1), data2(:,2)/100, config.n1_init, config.waveN_fit_min);

config.thk_init = mean([thk_init1, thk_init2]);

fprintf('文件1的FFT估算厚度: %.2f μm\n', thk_init1);
fprintf('文件2的FFT估算厚度: %.2f μm\n', thk_init2);
fprintf('平均厚度初值: %.2f μm\n', config.thk_init);


fprintf('\n');
disp('附件1（10°）：')
[result1.thk, result1.n2, result1.params] = process_file(data1, file1.angle, config);
fprintf('\n');
disp('附件1（15°）：')
[result2.thk, result2.n2, result2.params] = process_file(data2, file2.angle, config);


analyze_results(result1, result2, file1.angle, file2.angle);



%%
function [thk, n2, selParam] = process_file(data, incident_angle, config)
% 核心处理函数：数据预处理、拟合、结果输出


waveNum = data(:, 1);  % 全波段波数
R = data(:, 2) / 100;  % 全波段反射率
waveLen_full = 10000 ./ waveNum;

filter = waveNum > config.waveN_fit_min;
waveNum_fit = waveNum(filter);
R_fit = R(filter);
waveLen_fit = waveLen_full(filter);

x0 = [config.thk_init, config.n2_init, config.selParam_init];
lb = [config.thk_init*0.8, 2.0,  0.1, 0.001, 0.0001, 0.0001, 0.1,  0.001];
ub = [config.thk_init*1.2, 3.5,  20,  10,    5,      2,      150,  20];

[x_optimal, R_squared_fit] = global_fit_all_parameters(x0, lb, ub, waveNum_fit, R_fit, waveLen_fit, incident_angle);

thk = x_optimal(1);
n2 = x_optimal(2);
selParam = x_optimal(3:8);

fprintf('  - 最终厚度: %.2f μm\n', thk);
fprintf('  - 最终衬底折射率 n2: %.3f\n', n2);
fprintf('  - 在拟合区域的拟合优度 R²: %.4f\n', config.waveN_fit_min, R_squared_fit);

n1_all = calculate_sellmeier_n(waveLen_full, selParam);
plot_fitting_results(waveNum, R, n1_all, thk, n2, incident_angle, config.waveN_fit_min);
plot_refractive_index(waveLen_full, waveNum, n1_all);
end


function [x_optimal, R_squared] = global_fit_all_parameters(x0, lb, ub, waveNum_fit, R_fit, waveLen_fit, incident_angle)
% 全局优化函数
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


function thk = fft_thk_estimate(waveNum, R, n_avg, waveN_fit_min)
% FFT厚度估算
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
    thk_axis = (0:N/2-1) * 10000 / (2 * n_avg * N * dk);
    search_range = find(thk_axis > 5 & thk_axis < 200);

    if isempty(search_range), thk = 20;  return; end

    [~, max_idx_in_range] = max(fft_power(search_range));
    max_idx_rough = search_range(max_idx_in_range);

    correctNum = 3;
    DatePower1 = 0; DatePower2 = 0;
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

    thk = (f_corrected - 1) * 10000 / (2 * n_avg * N * dk);
end



function analyze_results(res1, res2, angle1, angle2)
% 最终结果对比与可靠性分析
    disp('最终结果对比：');
    fprintf('文件1 (%d°) -> 厚度: %.2f μm, 衬底折射率 n2: %.3f\n', angle1, res1.thk, res1.n2);
    fprintf('文件2 (%d°) -> 厚度: %.2f μm, 衬底折射率 n2: %.3f\n', angle2, res2.thk, res2.n2);

    thk1 = res1.thk; thk2 = res2.thk;
    params1 = res1.params; params2 = res2.params;

    thk_diff_percent = abs(thk1-thk2)/mean([thk1,thk2])*100;
    fprintf('两个角度测得的厚度分别为 %.2f μm 和 %.2f μm，相对差异为 %.2f%%。\n', thk1, thk2, thk_diff_percent);
    
    if thk_diff_percent < 5
        fprintf('  -> 结论: 一致性良好，结果可靠。\n');
        fprintf('  -> 推荐厚度值: %.2f ± %.2f μm\n', mean([thk1, thk2]), std([thk1, thk2]));
    else
        fprintf('  -> 结论: 差异较大。\n');
    end

    fprintf('\n 塞尔迈耶尔参数稳定性分析 \n');
    param_names = {'B1', 'B2', 'B3', 'C1', 'C2', 'C3'};
    is_stable = true;

    for i = 1:6
        diff_val = abs(params1(i) - params2(i))/mean([params1(i), params2(i)])*100;
        if diff_val > 50, is_stable = false; end
        fprintf('  %s: %.4f vs %.4f (差异: %.1f%%)\n', param_names{i}, params1(i), params2(i), diff_val);
    end

    if is_stable
        fprintf('  -> 结论: 参数在不同角度下拟合结果基本稳定。\n');
    else
        fprintf('  -> 结论: 部分参数差异较大，但只要最终的 n(λ) 曲线和厚度 d 一致。\n'); 
    end
end



%% 物理模型与绘图函数 

function R = model_R(params, ~, waveLen, theta0) 
    thk = params(1);
    n2 = params(2);
    selParam = params(3:8);
    n1 = calculate_sellmeier_n(waveLen, selParam);
    waveNum_model = 10000 ./ waveLen;
    R = model_R_vectorized([thk, n2], waveNum_model, n1, theta0);
end

function R = model_R_vectorized(params, waveNum, n1, theta0)
    thk = params(1); 
    n2 = params(2); 
    n0 = 1.0;
    n1 = n1(:); 
    waveNum = waveNum(:);

    sin_theta1 = n0 * sin(theta0) ./ n1; cos_theta1 = real(sqrt(1 - sin_theta1.^2));
    sin_theta2 = n0 * sin(theta0) / n2; cos_theta2 = real(sqrt(1 - sin_theta2^2));

    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    r12_s = (n1.*cos_theta1 - n2*cos_theta2) ./ (n1.*cos_theta1 + n2*cos_theta2);
    r01_p = (n1*cos(theta0) - n0.*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    r12_p = (n2*cos_theta1 - n1.*cos_theta2) ./ (n2*cos_theta1 + n1.*cos_theta2);

    delta = 4 * pi * n1 .* thk .* cos_theta1 .* waveNum / 10000;

    R_s = abs((r01_s + r12_s .* exp(1i*delta)) ./ (1 + r01_s .* r12_s .* exp(1i*delta))).^2;
    R_p = abs((r01_p + r12_p .* exp(1i*delta)) ./ (1 + r01_p .* r12_p .* exp(1i*delta))).^2;

    R = (R_s + R_p) / 2;
    R = real(R(:));
end

function n = calculate_sellmeier_n(waveLen, params)
    B1 = params(1); 
    B2 = params(2); 
    B3 = params(3);
    C1 = params(4); 
    C2 = params(5); 
    C3 = params(6);

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
    figure('Name', '拟合得到的折射率色散曲线', 'Position', [100, 100, 1200, 500]);

    subplot(1,2,1); plot(waveLen, n1_all, 'b-', 'LineWidth', 2); xlabel('波长 (μm)'); ylabel('外延层折射率 n_1'); title('折射率 vs 波长'); grid on;
    subplot(1,2,2); plot(waveNum, n1_all, 'r-', 'LineWidth', 2); xlabel('波数 (cm^{-1})'); ylabel('外延层折射率 n_1'); title('折射率 vs 波数'); grid on; xlim([min(waveNum), max(waveNum)]);
end


function plot_fitting_results(waveNum, R, n1, thk, n2, incident_angle, waveN_fit_min)
    theta0 = incident_angle * pi / 180;
    R_fitted = model_R_vectorized([thk, n2], waveNum, n1, theta0);
    R_squared_full = 1 - sum((R - R_fitted).^2) / sum((R - mean(R)).^2);

    figure('Name', ['拟合结果分析 (入射角 ' num2str(incident_angle) '°)'], 'Position', [150, 150, 1200, 600]);
    plot(waveNum, R*100, 'b.', 'MarkerSize', 5, 'DisplayName', '实测数据 ');
    hold on;
    plot(waveNum, R_fitted*100, 'r-', 'LineWidth', 2, 'DisplayName', '拟合模型');
    ylim_vals = ylim;

    h = fill([waveN_fit_min, max(waveNum), max(waveNum), waveN_fit_min], [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], 'k', 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'DisplayName', '拟合与处理区域');
    uistack(h, 'bottom');
    xlabel('波数 (cm^{-1})'); ylabel('反射率 (%)');
    title(sprintf('反射率拟合结果 (θ=%d°, d=%.2f μm, n_2=%.3f, 全局R²=%.4f)', incident_angle, thk, n2, R_squared_full));
    legend('Location', 'best'); grid on; xlim([min(waveNum), max(waveNum)]);
end