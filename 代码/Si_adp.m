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

disp('一、FFT计算厚度初值');
data1 = readmatrix(file1.path);
thk_init1 = fft_thk_estimate(data1(:,1), data1(:,2)/100, config.n1_init, 2000, file1.angle);
data2 = readmatrix(file2.path);
thk_init2 = fft_thk_estimate(data2(:,1), data2(:,2)/100, config.n1_init, 2000, file2.angle);
config.thk_init = mean([thk_init1, thk_init2]);
fprintf('10: %.2f μm\n', thk_init1);
fprintf('15: %.2f μm\n', thk_init2);
fprintf('平均厚度初值: %.2f μm\n', config.thk_init);

fprintf('\n二、双光束干涉精细拟合\n')
fprintf('\n10:\n')
[res1] = process(data1, file1.angle, config, file1.path);
fprintf('\n15:\n')
[res2] = process(data2, file2.angle, config, file2.path);

analyze_res(res1, res2, file1.angle, file2.angle);


function [res] = process(data, angle, config, output_filepath)
    waveNum = data(:, 1);
    R = data(:, 2) / 100;
    waveLen_full = 10000 ./ waveNum;
    filter = waveNum > config.waveN_fit_min;
    R_fit = R(filter);
    waveLen_fit = waveLen_full(filter);
    
    %% 柯西模型拟合
    disp('柯西模型拟合');
    x0_cauchy = [config.thk_init, config.cauchyParam_init, config.k1Param_A_init, config.k1Param_B_init, config.drudeParam_init];
    lb_cauchy = [config.thk_init*0.8, 3.35, 0, 0, 0, 100, 10];
    ub_cauchy = [config.thk_init*1.2, 3.6, 0.2, 1e-3, 4, 4000, 1000];
    
    [x_optimal_cauchy, R_squared_cauchy] = global_fit(x0_cauchy, lb_cauchy, ub_cauchy, R_fit, waveLen_fit, angle, config.epsilon_inf, 'cauchy');
    
    cauchy_res.thk = x_optimal_cauchy(1);
    cauchy_res.cauchyParam = x_optimal_cauchy(2:3);
    cauchy_res.k1Param = x_optimal_cauchy(4:5);
    cauchy_res.drudeParam = x_optimal_cauchy(6:7);
    cauchy_res.R_squared = R_squared_cauchy;
    cauchy_res.model_type = 'Cauchy';
    
    cauchy_res.n1_complex = cal_n_cauchy(waveLen_full, cauchy_res.cauchyParam, cauchy_res.k1Param);
    cauchy_res.n2_complex_full = cal_n_drude(waveNum, cauchy_res.drudeParam, config.epsilon_inf);
    
    fprintf('   R^2: %.6f\n', R_squared_cauchy);
    fprintf('   厚度: %.2f μm\n', cauchy_res.thk);
    fprintf('   柯西参数: A=%.3f, B=%.3f\n', cauchy_res.cauchyParam(1), cauchy_res.cauchyParam(2));
    fprintf('   k1参数: A=%.2e, B=%.2f\n', cauchy_res.k1Param(1), cauchy_res.k1Param(2));
    fprintf('   n2:Drude参数: ν_p=%.1f, Γ=%.1f\n', cauchy_res.drudeParam(1), cauchy_res.drudeParam(2));
    
    %% sel模型拟合
    disp('sel模型拟合');
    
    thk_cauchy = x_optimal_cauchy(1);
    k1_cauchy = x_optimal_cauchy(4:5);
    drude_cauchy = x_optimal_cauchy(6:7);
    
    x0_sel = [thk_cauchy, config.selParam_init, k1_cauchy, drude_cauchy];
    
    thk_rf = 0.03;
    lb_thk = thk_cauchy * (1 - thk_rf);
    ub_thk = thk_cauchy * (1 + thk_rf);
    
    sel_rf = 0.2;
    lb_sel_params = config.selParam_init * (1 - sel_rf);
    ub_sel_params = config.selParam_init * (1 + sel_rf);
    
    k1_rf = 0.3;
    lb_k1 = k1_cauchy * (1 - k1_rf);
    ub_k1 = k1_cauchy * (1 + k1_rf);
    
    drude_rf = 0.1;
    lb_drude = drude_cauchy * (1 - drude_rf);
    ub_drude = drude_cauchy * (1 + drude_rf);
    
    lb_sel = [lb_thk, lb_sel_params, lb_k1, lb_drude];
    ub_sel = [ub_thk, ub_sel_params, ub_k1, ub_drude];
    
    lb_sel = max(lb_sel, [config.thk_init*0.7, 0.1, 0.001, 0.0001, 0.0001, 0.1, 1, 0, 0, 50, 5]);
    ub_sel = min(ub_sel, [config.thk_init*1.3, 20, 5, 10, 2, 15, 2000, 1e-3, 4, 5000, 1500]);
    
    if lb_sel(8) < 0
        lb_sel(8) = 0;
    end
    if lb_sel(9) < 0
        lb_sel(9) = 0;
    end
    
    [x_optimal_sel, R_squared_sel] = global_fit(x0_sel, lb_sel, ub_sel, R_fit, waveLen_fit, angle, config.epsilon_inf, 'sellmeier');

    sel_res.thk = x_optimal_sel(1);
    sel_res.selParam = x_optimal_sel(2:7);
    sel_res.k1Param = x_optimal_sel(8:9);
    sel_res.drudeParam = x_optimal_sel(10:11);
    sel_res.R_squared = R_squared_sel;
    sel_res.model_type = 'sel';
    
    sel_res.n1_complex = cal_n_sellmeier(waveLen_full, sel_res.selParam, sel_res.k1Param);
    sel_res.n2_complex_full = cal_n_drude(waveNum, sel_res.drudeParam, config.epsilon_inf);
    
   ...
    fprintf('   R^2: %.6f\n', R_squared_sel);
    fprintf('   厚度: %.2f μm\n', sel_res.thk);
    fprintf('   赛尔迈耶参数: %.4f, %.4f, %.4f, %.4f, %.4f, %.1f\n', sel_res.selParam(1), sel_res.selParam(2), sel_res.selParam(3), sel_res.selParam(4), sel_res.selParam(5), sel_res.selParam(6));
    fprintf('   n2 Drude 参数: ν_p= %.1f , Γ= %.1f \n', sel_res.drudeParam(1), sel_res.drudeParam(2));
    ... 
    %% 选择最优模型
    disp('模型比较与选择');
    fprintf('   柯西模型 R^2 = %.6f\n', R_squared_cauchy);
    fprintf('   sel模型 R^2 = %.6f\n', R_squared_sel);
    
    if R_squared_sel > R_squared_cauchy
        res = sel_res;
        res.selected_model = 'sel';
        fprintf('   选择 sel 模型 (R^2提升: %.6f)\n', R_squared_sel - R_squared_cauchy);
    else
        res = cauchy_res;
        res.selected_model = 'Cauchy';
        fprintf('   选择 柯西 模型 (R^2更高或相等)\n');
    end
    
    % 保存两个模型的结果以供后续分析
    res.cauchy_res = cauchy_res;
    res.sel_res = sel_res;
    
    % 输出最终选择的模型参数
    fprintf('\n最终选择: %s 模型\n', res.selected_model);
    fprintf('   最终厚度: %.2f μm\n', res.thk);
    fprintf('   最终 R^2 = %.6f\n', res.R_squared);
    
    if strcmp(res.selected_model, 'sel')
        fprintf('   外延层 n1 (在 6μm): %.3f + %.4fi\n', ...
            real(interp1(waveLen_full, res.n1_complex, 6)), ...
            imag(interp1(waveLen_full, res.n1_complex, 6)));
    else
        fprintf('   外延层柯西参数: A=%.3f, B=%.3f\n', res.cauchyParam(1), res.cauchyParam(2));
        fprintf('   外延层 n1 (在 6μm): %.3f + %.4fi\n', ...
            real(interp1(waveLen_full, res.n1_complex, 6)), ...
            imag(interp1(waveLen_full, res.n1_complex, 6)));
    end
    fprintf('   k1参数: A=%.2e, B=%.2f\n', res.k1Param(1), res.k1Param(2));
    fprintf('   Drude参数: ν_p=%.1f , Γ=%.1f \n', res.drudeParam(1), res.drudeParam(2));
    
    % 保存结果到Excel
    try
        theta0_rad = angle * pi / 180;
        R_fit_full = compute_R(res.thk, res.n1_complex, res.n2_complex_full, waveNum, theta0_rad, 2);
        R_fit_per = R_fit_full * 100;
        writematrix(R_fit_per, output_filepath, 'Sheet', 1, 'Range', 'C2');
        fprintf('   成功将拟合结果保存到: %s (使用%s模型)\n', output_filepath, res.selected_model);
    catch ME
        fprintf('   保存Excel文件时出错: %s\n', ME.message);
    end
    
    figure;

    % 柯西
    subplot(1,2,1);
    R_fit_cauchy = compute_R(cauchy_res.thk, cauchy_res.n1_complex, cauchy_res.n2_complex_full, waveNum, angle * pi / 180, 2);
    plot(waveNum(2:end), R(2:end)*100, '-', 'Color', '#d74f44', 'LineWidth', 1.5, 'DisplayName', '实验数据');
    hold on;
    plot(waveNum, R_fit_cauchy*100, '-', 'Color', '#008ede', 'LineWidth', 1.3, 'DisplayName', sprintf('柯西模型 (R^2=%.4f)', cauchy_res.R_squared));
    xlabel('波数 (cm^{-1})');
    ylabel('反射率 (%)');
    legend('Location', 'best');
    grid on;
    
    % sel
    subplot(1,2,2);
    R_fit_sel = compute_R(sel_res.thk, sel_res.n1_complex, sel_res.n2_complex_full, waveNum, angle * pi / 180, 2);
    plot(waveNum(2:end), R(2:end)*100, '-', 'Color', '#d74f44', 'LineWidth', 1.5, 'DisplayName', '实验数据');
    hold on;
    plot(waveNum, R_fit_sel*100, '-', 'Color', '#3f8819', 'LineWidth', 1.3, 'DisplayName', sprintf('sel模型 (R^2=%.4f)', sel_res.R_squared));
    xlabel('波数 (cm^{-1})');
    ylabel('反射率 (%)');
    legend('Location', 'best');
    grid on;
        
    % 绘制最优模型的光学常数
    plot_optical_constants(waveLen_full, waveNum, res.n1_complex, sprintf('外延层 (%s)', res.selected_model));
    plot_optical_constants(waveLen_full, waveNum, res.n2_complex_full, '衬底 (Drude)');
end

function analyze_res(res1, res2, angle1, angle2)
    disp('最终结果对比');
    fprintf('文件3 (%d°) -> 厚度: %.2f μm, ν_p: %.1f , Γ: %.1f\n ', angle1, res1.thk, res1.drudeParam(1), res1.drudeParam(2));
    fprintf('文件4 (%d°) -> 厚度: %.2f μm, ν_p: %.1f , Γ: %.1f\n ', angle2, res2.thk, res2.drudeParam(1), res2.drudeParam(2));
    thk1 = res1.thk; thk2 = res2.thk;
    thk_diff_percent = abs(thk1-thk2)/mean([thk1,thk2])*100;
    fprintf('两个角度测得的厚度分别为 %.2f μm 和 %.2f μm，相对差异为 %.2f%%。\n', thk1, thk2, thk_diff_percent);
    if thk_diff_percent < 10
        fprintf('   -> 结论: 厚度一致性良好，结果可靠。\n');
        fprintf('   -> 推荐厚度值: %.2f ± %.2f μm\n', mean([thk1, thk2]), std([thk1, thk2]));
    else
        fprintf('   -> 结论: 厚度差异较大。\n');
    end
end
