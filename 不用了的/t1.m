clear; % 清除工作区变量
clc; % 清除命令行窗口
close all; % 关闭所有图形窗口

filename = '附件2.xlsx';
 dataTable = readtable(filename);

x = dataTable{:, 1}; % 提取第一列: 波数 (cm^-1)
y = dataTable{:, 2}; % 提取第二列: 反射率 (%)

figure;
plot(x, y, 'b-', 'LineWidth', 1.5);
title('原始数据: 反射率 vs 波数');
xlabel('波数 (cm^{-1})');
ylabel('反射率 (%)');
grid on;
set(gca, 'FontSize', 12);

% idx = (x>=400)&(x<730);
idx = x >= 1000; % 找到下标 ≥ 1000 cm^-1 的位置
x_cut = x(idx); % 截取后的波数
y_cut = y(idx); % 截取后的反射率

win = hann(length(y_cut));
y_win = y_cut .* win; 

dx = mean(diff(x_cut)); % 平均采样间隔(cm^-1)
Fs = 1 / dx; % "采样频率"
N = length(y_win); % 信号长度
Y = fft(y_win); % 加窗后的 FFT
P2 = abs(Y / N); % 双边振幅
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2 * P1(2:end-1);
% 空间域坐标（单位 cm）
f = Fs * (0:(N/2)) / N;


% 找到最大峰值
threshold_f = 0.001;           % 频率阈值 (cm)
valid_idx = find(f > threshold_f); % 找到所有 f > 0.01 cm 的索引

% --- 在 f > 0.01 cm 范围内找最大值 ---
[maxVal, peakIdx] = max(P1(valid_idx)); % 注意：peakIdx 是 valid_idx 的局部索引
peakIdx = valid_idx(peakIdx); % 转换为全局索引
f_peak = f(peakIdx); % 原始峰值频率

% 能量重心校正
correctNum = 5; % 校正范围 (可调)
[VPP_corrected, f_corrected] = ADC_FFT_Get_Wave_Mes(peakIdx, Fs, P1, correctNum, N);

% 显示结果
disp(['原始主峰频率 = ', num2str(f_peak), ' cm']);
disp(['能量中心校正频率 = ', num2str(f_corrected), ' cm']);
% disp(['校正后峰峰值 = ', num2str(VPP_corrected)]);


figure;
plot(f, P1, 'r-', 'LineWidth', 1.5);
hold on;
xlabel('空间位置 (cm)');
ylabel('振幅');

grid on;
set(gca, 'FontSize', 12);
xlim([0, max(f)]);
