%% 碳化硅外延层厚度测定主程序
clear; close all; clc;

%% 1. 读取数据
% 读取附件1和附件2的数据（10度和15度入射角）
data1 = readmatrix('附件1.xlsx'); % 10度入射角
data2 = readmatrix('附件2.xlsx'); % 15度入射角

% 提取波数和反射率
wavenumber1 = data1(:,1); % 波数 (cm^-1)
reflectance1 = data1(:,2); % 反射率 (%)
wavenumber2 = data2(:,1);
reflectance2 = data2(:,2);

%% 2. 数据预处理 - 截取波数大于1000的数据
idx1 = wavenumber1 >= 1000;
k1 = wavenumber1(idx1);
R1 = reflectance1(idx1);

idx2 = wavenumber2 >= 1000;
k2 = wavenumber2(idx2);
R2 = reflectance2(idx2);

%% 3. 使用塞尔迈耶尔方程拟合外延层折射率n1
% 塞尔迈耶尔方程: n^2 = 1 + Σ(Ai*λ^2)/(λ^2 - Bi)
% 其中 λ = 10000/k (转换为微米)

% 选择波数大于2000的数据进行拟合
idx_fit = k1 > 2000;
k_fit = k1(idx_fit);
lambda_fit = 10000 ./ k_fit; % 转换为波长(微米)

% 碳化硅折射率的初始参考值（根据文献）
n_ref = 2.65; % SiC典型折射率

% 构建塞尔迈耶尔方程的拟合函数
sellmeier = @(params, lambda) sqrt(1 + ...
    params(1)*lambda.^2./(lambda.^2 - params(4)) + ...
    params(2)*lambda.^2./(lambda.^2 - params(5)) + ...
    params(3)*lambda.^2./(lambda.^2 - params(6)));

% 初始参数估计 [A1, A2, A3, B1, B2, B3]
initial_params = [5.0, 0.5, 0.1, 0.1, 0.5, 10];

% 使用最小二乘法拟合
% 需要从反射率数据中提取折射率信息
% 对于单层反射，反射率 R = ((n1-n0)/(n1+n0))^2，其中n0=1(空气)
% 因此 n1 = (1+sqrt(R/100))/(1-sqrt(R/100))

% 估算高频区域的折射率（干涉效应较小）
R_high = R1(idx_fit);
n_estimated = (1 + sqrt(R_high/100)) ./ (1 - sqrt(R_high/100));

% 拟合塞尔迈耶尔参数
options = optimset('Display','iter','TolFun',1e-8,'TolX',1e-8);
params_fit = lsqcurvefit(sellmeier, initial_params, lambda_fit, n_estimated, ...
    [0.1, 0.01, 0.001, 0.01, 0.1, 1], ...  % 下界
    [10, 5, 1, 1, 5, 50], ...               % 上界
    options);

% 提取拟合参数
A1 = params_fit(1); A2 = params_fit(2); A3 = params_fit(3);
B1 = params_fit(4); B2 = params_fit(5); B3 = params_fit(6);

fprintf('塞尔迈耶尔参数拟合结果:\n');
fprintf('A1 = %.6f, A2 = %.6f, A3 = %.6f\n', A1, A2, A3);
fprintf('B1 = %.6f, B2 = %.6f, B3 = %.6f\n', B1, B2, B3);

% 计算所有波数对应的折射率n1
lambda_all = 10000 ./ k1; % 所有波长
n1_all = sellmeier(params_fit, lambda_all);

%% 4. FFT分析找到主周期
% 加窗处理
win = hann(length(R1));
R1_windowed = R1 .* win;

% 计算采样间隔
dk = mean(diff(k1)); % 平均波数间隔
Fs = 1 / dk;         % 采样频率

% FFT分析
N = length(R1_windowed);
Y = fft(R1_windowed);
P2 = abs(Y / N);
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2 * P1(2:end-1);

% 频率坐标（单位：cm）
f = Fs * (0:(N/2)) / N;

% 找到最大峰值
threshold_f = 0.001; % 频率阈值（cm）
valid_idx = find(f > threshold_f);
[maxVal, peakIdx_temp] = max(P1(valid_idx));
peakIdx = valid_idx(peakIdx_temp);
f_peak = f(peakIdx);

% 能量重心校正（更精确的峰值位置）
correctNum = 5; % 校正范围
[f_corrected, ~] = energy_centroid_correction(peakIdx, f, P1, correctNum);

fprintf('\n原始主频率 = %.4f cm\n', f_peak);
fprintf('能量重心校正后频率 = %.4f cm\n', f_corrected);

%% 5. 衬底折射率n2的估计
% 使用中心频率处的折射率作为衬底折射率的初始估计
k_center = mean(k1);
lambda_center = 10000 / k_center;
n2_init = sellmeier(params_fit, lambda_center) * 0.95; % 衬底折射率通常略小

%% 6. 外延层厚度的初始估计
% 对于垂直入射，光程差 = 2*n1*d
% 干涉条件：2*n1*d = m*λ，其中m是干涉级次
% 在波数域：2*n1*d*k = m*10000
% 周期 Δk = 10000/(2*n1*d)
% 因此 d = 10000/(2*n1*Δk)

% 使用平均折射率
n1_avg = mean(n1_all);
d_init = 1 / (2 * n1_avg * f_corrected); % 初始厚度估计（微米）

fprintf('\n初始厚度估计 d = %.2f μm\n', d_init * 10000);

%% 7. 精确拟合：考虑入射角和完整的干涉模型
theta1 = 10 * pi/180; % 入射角（弧度）
theta2 = 15 * pi/180;

% 干涉反射率模型（考虑入射角）
% R = |r1 + r2*exp(2i*beta)|^2 / |1 + r1*r2*exp(2i*beta)|^2
% 其中 beta = 2*pi*n1*d*cos(theta_r)*k/10000
% theta_r是折射角，由Snell定律确定

reflectance_model = @(params, k, theta) calculate_reflectance(params, k, theta, params_fit);

% 参数：[d(微米), n2, Δn(外延层与衬底折射率差的修正)]
initial_params2 = [d_init*10000, n2_init, 0];

% 同时拟合两个角度的数据
k_combined = [k1; k2];
R_combined = [R1; R2];
theta_combined = [ones(size(k1))*theta1; ones(size(k2))*theta2];

% 非线性最小二乘拟合
options2 = optimoptions('lsqcurvefit','Display','iter','MaxFunctionEvaluations',5000);
params_final = lsqcurvefit(@(p,k) reflectance_model(p, k_combined, theta_combined), ...
    initial_params2, k_combined, R_combined, ...
    [d_init*10000*0.5, n2_init*0.9, -0.1], ...  % 下界
    [d_init*10000*1.5, n2_init*1.1, 0.1], ...   % 上界
    options2);

d_final = params_final(1);      % 厚度（微米）
n2_final = params_final(2);     % 衬底折射率
delta_n = params_final(3);      % 折射率修正

fprintf('\n最终拟合结果:\n');
fprintf('外延层厚度 d = %.3f μm\n', d_final);
fprintf('衬底折射率 n2 = %.4f\n', n2_final);
fprintf('折射率修正 Δn = %.4f\n', delta_n);

%% 8. 结果可视化
figure('Position', [100, 100, 1200, 800]);

% 子图1：实测与拟合对比（10度）
subplot(2,2,1);
R_fit1 = reflectance_model(params_final, k1, theta1*ones(size(k1)));
plot(k1, R1, 'b-', 'LineWidth', 1, 'DisplayName', '实测数据');
hold on;
plot(k1, R_fit1, 'r--', 'LineWidth', 1, 'DisplayName', '拟合结果');
xlabel('波数 (cm^{-1})');
ylabel('反射率 (%)');
title('10°入射角');
legend('Location', 'best');
grid on;

% 子图2：实测与拟合对比（15度）
subplot(2,2,2);
R_fit2 = reflectance_model(params_final, k2, theta2*ones(size(k2)));
plot(k2, R2, 'b-', 'LineWidth', 1, 'DisplayName', '实测数据');
hold on;
plot(k2, R_fit2, 'r--', 'LineWidth', 1, 'DisplayName', '拟合结果');
xlabel('波数 (cm^{-1})');
ylabel('反射率 (%)');
title('15°入射角');
legend('Location', 'best');
grid on;

% 子图3：FFT频谱
subplot(2,2,3);
plot(f, P1, 'b-', 'LineWidth', 1);
hold on;
plot(f_peak, P1(peakIdx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
plot(f_corrected, P1(peakIdx), 'g^', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('空间频率 (cm)');
ylabel('幅度');
title('FFT频谱分析');
legend('FFT频谱', '原始峰值', '校正后峰值', 'Location', 'best');
xlim([0, 0.01]);
grid on;

% 子图4：折射率随波长变化
subplot(2,2,4);
plot(lambda_all, n1_all, 'b-', 'LineWidth', 2);
hold on;
yline(n2_final, 'r--', 'LineWidth', 1.5);
xlabel('波长 (μm)');
ylabel('折射率');
title('折射率分布');
legend('外延层 n_1', '衬底 n_2', 'Location', 'best');
grid on;

%% 9. 误差分析
% 计算拟合残差
residual1 = R1 - R_fit1;
residual2 = R2 - R_fit2;
RMSE1 = sqrt(mean(residual1.^2));
RMSE2 = sqrt(mean(residual2.^2));

fprintf('\n拟合误差分析:\n');
fprintf('10°入射角 RMSE = %.4f%%\n', RMSE1);
fprintf('15°入射角 RMSE = %.4f%%\n', RMSE2);

%% 辅助函数
function R = calculate_reflectance(params, k, theta, sellmeier_params)
    % 计算反射率
    d = params(1);      % 厚度（微米）
    n2 = params(2);     % 衬底折射率
    delta_n = params(3); % 折射率修正
    
    % 计算外延层折射率
    lambda = 10000 ./ k;
    n1 = sqrt(1 + ...
        sellmeier_params(1)*lambda.^2./(lambda.^2 - sellmeier_params(4)) + ...
        sellmeier_params(2)*lambda.^2./(lambda.^2 - sellmeier_params(5)) + ...
        sellmeier_params(3)*lambda.^2./(lambda.^2 - sellmeier_params(6))) + delta_n;
    
    % 空气折射率
    n0 = 1;
    
    % Snell定律计算折射角
    sin_theta_r = sin(theta) .* n0 ./ n1;
    cos_theta_r = sqrt(1 - sin_theta_r.^2);
    
    % 菲涅尔反射系数（s偏振）
    r01 = (n0*cos(theta) - n1.*cos_theta_r) ./ (n0*cos(theta) + n1.*cos_theta_r);
    r12 = (n1.*cos_theta_r - n2.*cos_theta_r) ./ (n1.*cos_theta_r + n2.*cos_theta_r);
    
    % 相位
    beta = 2*pi*n1.*d.*cos_theta_r.*k/10000;
    
    % 总反射率
    numerator = abs(r01 + r12.*exp(2i*beta)).^2;
    denominator = abs(1 + r01.*r12.*exp(2i*beta)).^2;
    R = 100 * numerator ./ denominator;
end

function [f_corrected, P_corrected] = energy_centroid_correction(peakIdx, f, P, correctNum)
    % 能量重心校正
    idx_range = max(1, peakIdx-correctNum):min(length(f), peakIdx+correctNum);
    f_range = f(idx_range);
    P_range = P(idx_range);
    
    % 计算能量重心
    f_corrected = sum(f_range .* P_range) / sum(P_range);
    P_corrected = interp1(f, P, f_corrected);
end