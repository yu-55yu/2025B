clc;
clear;
close all;

file1.path = '附件3.xlsx';
file1.angle = 10;
file2.path = '附件4.xlsx';
file2.angle = 15;

config.waveN_fit_min = 400;
config.n1_init = 3.48;
config.cauchyParam_init = [3.42, 0.05];
config.selParam_init = [10.6684293, 0.0030434748, 1.5413340, 0.301516485, 1.13475115, 1104];
config.k1Param_A_init = 1e-5;
config.k1Param_B_init = 2.0;
config.drudeParam_init = [2000, 200];
config.epsilon_inf = 11.7;

disp('计算厚度初值');
data1 = readmatrix(file1.path);
thk_init1 = fft_thk_estimate(data1(:,1), data1(:,2)/100, config.n1_init, 2000, file1.angle);
data2 = readmatrix(file2.path);
thk_init2 = fft_thk_estimate(data2(:,1), data2(:,2)/100, config.n1_init, 2000, file2.angle);
config.thk_init = mean([thk_init1, thk_init2]);
fprintf('文件3 (10°) FFT估算厚度: %.2f μm\n', thk_init1);
fprintf('文件4 (15°) FFT估算厚度: %.2f μm\n', thk_init2);
fprintf('平均厚度初值: %.2f μm\n', config.thk_init);

disp('--- 开始处理附件3 (10°) ---')
[result1] = process(data1, file1.angle, config, file1.path);
disp('--- 开始处理附件4 (15°) ---')
[result2] = process(data2, file2.angle, config, file2.path);

analyze_res(result1, result2, file1.angle, file2.angle);


function [result] = process(data, angle, config, output_filepath)
    waveNum = data(:, 1);
    R = data(:, 2) / 100;
    waveLen_full = 10000 ./ waveNum;
    filter = waveNum > config.waveN_fit_min;
    R_fit = R(filter);
    waveLen_fit = waveLen_full(filter);
    
    %% Phase 1: 柯西模型拟合
    disp('--- Phase 1: 柯西模型拟合');
    x0_cauchy = [config.thk_init, config.cauchyParam_init, config.k1Param_A_init, config.k1Param_B_init, config.drudeParam_init];
    lb_cauchy = [config.thk_init*0.8, 3.35, 0, 0, 0, 100, 10];
    ub_cauchy = [config.thk_init*1.2, 3.6, 0.2, 1e-3, 4, 4000, 1000];
    
    [x_optimal_cauchy, R_squared_cauchy] = global_fit(x0_cauchy, lb_cauchy, ub_cauchy, R_fit, waveLen_fit, angle, config.epsilon_inf, 'cauchy');
    
    % 保存柯西模型结果
    cauchy_result.thk = x_optimal_cauchy(1);
    cauchy_result.cauchyParam = x_optimal_cauchy(2:3);
    cauchy_result.k1Param = x_optimal_cauchy(4:5);
    cauchy_result.drudeParam = x_optimal_cauchy(6:7);
    cauchy_result.R_squared = R_squared_cauchy;
    cauchy_result.model_type = 'Cauchy';
    
    % 计算柯西模型的复折射率
    cauchy_result.n1_complex_full = cal_n_cauchy(waveLen_full, cauchy_result.cauchyParam, cauchy_result.k1Param);
    cauchy_result.n2_complex_full = cal_n_drude(waveNum, cauchy_result.drudeParam, config.epsilon_inf);
    
    fprintf('   柯西模型完成. R^2 = %.6f\n', R_squared_cauchy);
    fprintf('   柯西拟合厚度: %.2f μm\n', cauchy_result.thk);
    fprintf('   柯西参数: A=%.3f, B=%.3f\n', cauchy_result.cauchyParam(1), cauchy_result.cauchyParam(2));
    fprintf('   k1参数: A=%.2e, B=%.2f\n', cauchy_result.k1Param(1), cauchy_result.k1Param(2));
    fprintf('   Drude参数: ν_p=%.1f, Γ=%.1f\n', cauchy_result.drudeParam(1), cauchy_result.drudeParam(2));
    
    %% Phase 2: Sellmeier模型拟合
    disp('--- Phase 2: Sellmeier模型拟合');
    
    % 基于柯西拟合结果设置Sellmeier初始值
    thk_cauchy = x_optimal_cauchy(1);
    k1_cauchy = x_optimal_cauchy(4:5);
    drude_cauchy = x_optimal_cauchy(6:7);
    
    x0_sel = [thk_cauchy, config.selParam_init, k1_cauchy, drude_cauchy];
    
    % 基于柯西拟合结果设置更紧的边界
    thk_range_factor = 0.03;
    lb_thk = thk_cauchy * (1 - thk_range_factor);
    ub_thk = thk_cauchy * (1 + thk_range_factor);
    
    sel_range_factor = 0.2;
    lb_sel_params = config.selParam_init * (1 - sel_range_factor);
    ub_sel_params = config.selParam_init * (1 + sel_range_factor);
    
    k1_range_factor = 0.3;
    lb_k1 = k1_cauchy * (1 - k1_range_factor);
    ub_k1 = k1_cauchy * (1 + k1_range_factor);
    
    drude_range_factor = 0.1;
    lb_drude = drude_cauchy * (1 - drude_range_factor);
    ub_drude = drude_cauchy * (1 + drude_range_factor);
    
    lb_sel = [lb_thk, lb_sel_params, lb_k1, lb_drude];
    ub_sel = [ub_thk, ub_sel_params, ub_k1, ub_drude];
    
    % 确保边界在合理的物理范围内
    lb_sel = max(lb_sel, [config.thk_init*0.7, 0.1, 0.001, 0.0001, 0.0001, 0.1, 1, 0, 0, 50, 5]);
    ub_sel = min(ub_sel, [config.thk_init*1.3, 20, 5, 10, 2, 15, 2000, 1e-3, 4, 5000, 1500]);
    
    if lb_sel(8) < 0
        lb_sel(8) = 0;
    end
    if lb_sel(9) < 0
        lb_sel(9) = 0;
    end
    
    [x_optimal_sel, R_squared_sel] = global_fit(x0_sel, lb_sel, ub_sel, R_fit, waveLen_fit, angle, config.epsilon_inf, 'sellmeier');

    % 保存Sellmeier模型结果
    sellmeier_result.thk = x_optimal_sel(1);
    sellmeier_result.selParam = x_optimal_sel(2:7);
    sellmeier_result.k1Param = x_optimal_sel(8:9);
    sellmeier_result.drudeParam = x_optimal_sel(10:11);
    sellmeier_result.R_squared = R_squared_sel;
    sellmeier_result.model_type = 'Sellmeier';
    
    % 计算Sellmeier模型的复折射率
    sellmeier_result.n1_complex_full = cal_n_sellmeier(waveLen_full, sellmeier_result.selParam, sellmeier_result.k1Param);
    sellmeier_result.n2_complex_full = cal_n_drude(waveNum, sellmeier_result.drudeParam, config.epsilon_inf);
    
    fprintf('   Sellmeier模型完成. R^2 = %.6f\n', R_squared_sel);
    fprintf('   Sellmeier拟合厚度: %.2f μm\n', sellmeier_result.thk);
    fprintf('   最终衬底 Drude 参数 ν_p: %.1f cm??, Γ: %.1f cm??\n', sellmeier_result.drudeParam(1), sellmeier_result.drudeParam(2));
    
    %% 选择最优模型
    disp('--- 模型比较与选择 ---');
    fprintf('   柯西模型 R^2 = %.6f\n', R_squared_cauchy);
    fprintf('   Sellmeier模型 R^2 = %.6f\n', R_squared_sel);
    
    if R_squared_sel > R_squared_cauchy
        result = sellmeier_result;
        result.selected_model = 'Sellmeier';
        fprintf('   >>> 选择 Sellmeier 模型 (R^2提升: %.6f)\n', R_squared_sel - R_squared_cauchy);
    else
        result = cauchy_result;
        result.selected_model = 'Cauchy';
        fprintf('   >>> 选择 柯西 模型 (R^2更高或相等)\n');
    end
    
    % 保存两个模型的结果以供后续分析
    result.cauchy_result = cauchy_result;
    result.sellmeier_result = sellmeier_result;
    
    % 输出最终选择的模型参数
    fprintf('\n=== 最终选择: %s 模型 ===\n', result.selected_model);
    fprintf('   最终厚度: %.2f μm\n', result.thk);
    fprintf('   最终 R^2 = %.6f\n', result.R_squared);
    
    if strcmp(result.selected_model, 'Sellmeier')
        fprintf('   外延层 n1 (在 6μm): %.3f + %.4fi\n', ...
            real(interp1(waveLen_full, result.n1_complex_full, 6)), ...
            imag(interp1(waveLen_full, result.n1_complex_full, 6)));
    else
        fprintf('   外延层柯西参数: A=%.3f, B=%.3f\n', result.cauchyParam(1), result.cauchyParam(2));
        fprintf('   外延层 n1 (在 6μm): %.3f + %.4fi\n', ...
            real(interp1(waveLen_full, result.n1_complex_full, 6)), ...
            imag(interp1(waveLen_full, result.n1_complex_full, 6)));
    end
    fprintf('   k1参数: A=%.2e, B=%.2f\n', result.k1Param(1), result.k1Param(2));
    fprintf('   Drude参数: ν_p=%.1f cm??, Γ=%.1f cm??\n', result.drudeParam(1), result.drudeParam(2));
    
    % 保存结果到Excel
    try
        theta0_rad = angle * pi / 180;
        R_fit_full = compute_R(result.thk, result.n1_complex_full, result.n2_complex_full, waveNum, theta0_rad, 0);
        R_fit_per = R_fit_full * 100;
        writematrix(R_fit_per, output_filepath, 'Sheet', 1, 'Range', 'C2');
        fprintf('   成功将拟合结果保存到: %s (使用%s模型)\n', output_filepath, result.selected_model);
    catch ME
        fprintf('   保存Excel文件时出错: %s\n', ME.message);
    end
    
    % 绘图 - 显示两个模型的对比
    figure;

    % --- 柯西模型子图 ---
    subplot(1,2,1);
    R_fit_cauchy = compute_R(cauchy_result.thk, cauchy_result.n1_complex_full, cauchy_result.n2_complex_full, waveNum, angle * pi / 180, 0);

    plot(waveNum(2:end), R(2:end)*100, '-', 'Color', '#d74f44', 'LineWidth', 1.5, 'DisplayName', '实验数据');
    hold on;
    % 使用 #fc8a9e 颜色绘制柯西模型拟合数据
    plot(waveNum, R_fit_cauchy*100, '-', 'Color', '#008ede', 'LineWidth', 1.3, 'DisplayName', sprintf('柯西模型 (R^2=%.4f)', cauchy_result.R_squared));
    xlabel('波数 (cm^{-1})');
    ylabel('反射率 (%)');
    legend('Location', 'best');
    grid on;
    
    % --- Sellmeier模型子图 ---
    subplot(1,2,2);
    R_fit_sellmeier = compute_R(sellmeier_result.thk, sellmeier_result.n1_complex_full, sellmeier_result.n2_complex_full, waveNum, angle * pi / 180, 0);
    % 使用 #008ede 颜色绘制实验数据
    plot(waveNum(2:end), R(2:end)*100, '-', 'Color', '#d74f44', 'LineWidth', 1.5, 'DisplayName', '实验数据');
    hold on;
    % 使用 #3f8819 颜色绘制Sellmeier模型拟合数据
    plot(waveNum, R_fit_sellmeier*100, '-', 'Color', '#3f8819', 'LineWidth', 1.3, 'DisplayName', sprintf('Sellmeier模型 (R^2=%.4f)', sellmeier_result.R_squared));
    xlabel('波数 (cm^{-1})');
    ylabel('反射率 (%)');
    legend('Location', 'best');
    grid on;
        
    % 绘制最优模型的光学常数
    plot_optical_constants(waveLen_full, waveNum, result.n1_complex_full, sprintf('外延层 (%s)', result.selected_model));
    plot_optical_constants(waveLen_full, waveNum, result.n2_complex_full, '衬底 (Drude)');
end

function analyze_res(res1, res2, angle1, angle2)
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
