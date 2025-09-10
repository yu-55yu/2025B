clc;
clear;
close all;

file1.path = '附件1.xlsx';
file1.angle = 10; % 入射角 (度)
file2.path = '附件2.xlsx';
file2.angle = 15; % 入射角 (度)

config.waveN_fit_min = 1500; % 拟合起始波数 (cm-1)
config.n1_init = 2.58;       % 外延层折射率初值 (用于FFT估算厚度)
config.n2_real_init = 2.55;  % 衬底折射率初值
config.k1_A_init = 0.001;    % 外延层消光系数参数A初值
config.k1_B_init = 2.0;      % 外延层消光系数参数B初值

% sel模型 B1, B2, B3, C1, C2, C3 参数初值
config.selParam_init = [5.5, 0.2, 0.05, 0.027, 100, 0.01];



disp('FFT计算厚度初值');
data1 = readmatrix(file1.path);
thk_init1 = fft_thk_estimate(data1(:,1), data1(:,2)/100, config.n1_init, config.waveN_fit_min, file1.angle);
data2 = readmatrix(file2.path);
thk_init2 = fft_thk_estimate(data2(:,1), data2(:,2)/100, config.n1_init, config.waveN_fit_min, file2.angle);
config.thk_init = mean([thk_init1, thk_init2]);
fprintf('文件1 (10°) FFT估算厚度: %.2f μm\n', thk_init1);
fprintf('文件2 (15°) FFT估算厚度: %.2f μm\n', thk_init2);
fprintf('<strong>平均厚度初值: %.2f μm</strong>\n', config.thk_init);


disp(' ');
disp('10：');
[res1] = process(data1, file1.angle, config, file1.path);
disp('15：');
[res2] = process(data2, file2.angle, config, file2.path);


disp(' ');
disp('分析：');
analyze_res(res1, res2, file1.angle, file2.angle);




%% 数据处理和拟合主函数
function [res] = process(data, angle, config, output_filepath)
epsilon= 11.7;
waveNum = data(:, 1);       % 波数 (cm-1)
R = data(:, 2) / 100;       % 反射率 (0-1)
waveLen = 10000 ./ waveNum; % 波长 (μm)

filter = waveNum >= config.waveN_fit_min;
waveNum_fit = waveNum(filter);
waveLen_fit = 10000 ./ waveNum_fit;

R_fit = R(filter);

% 厚度, 衬底n2, sel(6), k参数(2)
x0 = [config.thk_init, config.n2_real_init, config.selParam_init, config.k1_A_init, config.k1_B_init];
lb = [config.thk_init*0.8, 2.0, 0.01, 0.001, 0.0001, 0.0001, 0.1, 0.001, 0, 0];
ub = [config.thk_init*1.2, 3.5, 20, 10, 5, 2, 150, 20, 0.01, 4];

model_type = 'real_substrate';
[x_optimal, R_squared_fit] = global_fit(x0, lb, ub,  R_fit,waveLen_fit, angle, epsilon ,model_type);

res.thk = x_optimal(1);
res.n2 = x_optimal(2);
res.selParam = x_optimal(3:8);
res.k1Param = x_optimal(9:10);
res.n1_complex = cal_n_sellmeier(waveLen, res.selParam, res.k1Param);

fprintf('    最终厚度: %.2f μm\n', res.thk);
fprintf('    最终外延层 n1 (在 6μm): %.3f + %.4fi\n', real(interp1(waveLen, res.n1_complex, 6)), imag(interp1(waveLen, res.n1_complex, 6)));
fprintf('    最终衬底 n2: %.3f\n', res.n2);
fprintf('    在拟合区域的拟合优度 R?: %.4f\n', R_squared_fit);

theta0_rad = angle * pi / 180;


R_fit_full = compute_R(res.thk, res.n1_complex, res.n2, waveNum, theta0_rad,2);


R_fit_per = R_fit_full * 100;
header = {'拟合反射率'};
writecell(header, output_filepath, 'Sheet', 1, 'Range', 'C1');
writematrix(R_fit_per, output_filepath, 'Sheet', 1, 'Range', 'C2');


plot_fit_res(waveNum, R, R_fit_full, angle, config.waveN_fit_min);
plot_optical_constants(waveLen, waveNum, res.n1_complex, '外延层 (sel)');
end


function analyze_res(res1, res2, angle1, angle2)
% ANALYZE_RES 分析和比较两个角度的拟合结果
%
% 输入:
%   res1 - 第一个角度的拟合结果结构体
%   res2 - 第二个角度的拟合结果结构体
%   angle1  - 第一个角度的值 (度)
%   angle2  - 第二个角度的值 (度)


fprintf('<strong> 关键物理参数对比 </strong>\n');

% 比较外延层厚度 (thk)
thk1 = res1.thk;
thk2 = res2.thk;
thk_diff_abs = abs(thk1 - thk2);
thk_diff_rel = thk_diff_abs / mean([thk1, thk2]) * 100;
fprintf('外延层厚度:\n');
fprintf('  %d° 拟合结果: %.3f μm\n', angle1, thk1);
fprintf('  %d° 拟合结果: %.3f μm\n', angle2, thk2);
fprintf('  -> 绝对差异: %.3f μm\n', thk_diff_abs);
fprintf('  -> 相对差异: %.2f %%\n', thk_diff_rel);
if thk_diff_rel < 10.0
    fprintf('  结论: 厚度结果具有高度一致性。\n\n');
else
    fprintf('  结论: 厚度结果存在一定差异，请检查。\n\n');
end

% 比较衬底折射率 (n2)
n2_1 = res1.n2;
n2_2 = res2.n2;
n2_diff_abs = abs(n2_1 - n2_2);
n2_diff_rel = n2_diff_abs / mean([n2_1, n2_2]) * 100;
fprintf('衬底折射率 (n2):\n');
fprintf('  %d° 拟合结果: %.3f\n', angle1, n2_1);
fprintf('  %d° 拟合结果: %.3f\n', angle2, n2_2);
fprintf('  -> 绝对差异: %.3f\n', n2_diff_abs);
fprintf('  -> 相对差异: %.2f %%\n', n2_diff_rel);
if n2_diff_rel < 10.0
    fprintf('  结论: 衬底折射率结果一致性良好。\n\n');
else
    fprintf('  结论: 衬底折射率结果存在一定差异，请检查。\n\n');
end
end