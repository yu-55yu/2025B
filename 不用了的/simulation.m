clc; clear; close all;

%% 参数设置
n0 = 1.0;          % 入射介质折射率（空气）
n1 = 3.5;          % 外延层折射率（假设常数）
n2 = 4.0;          % 衬底折射率
d  = 500e-9;       % 外延层厚度 (m)
theta0 = 0;        % 入射角 (弧度)
c = 3e8;

% 波数范围 (cm^-1)  转换成波长范围
nu_min = 1000;     % cm^-1
nu_max = 10000;     % cm^-1
N = 2^12;          % 采样点数
nu = linspace(nu_min,nu_max,N); % 波数 cm^-1
lambda = 1./(nu*100);           % 波长 (m)

% 在薄膜内的传播角度
theta1 = asin(n0/n1 * sin(theta0));

%% 情况1：两束干涉
% 幅值反射系数
r01 = (n0*cos(theta0)-n1*cos(theta1))/(n0*cos(theta0)+n1*cos(theta1));
r12 = (n1*cos(theta1)-n2*cos(theta0))/(n1*cos(theta1)+n2*cos(theta0));

% 只取两束相加：E = r01 + r12*exp(iδ)
delta = 4*pi*n1*d*cos(theta1)./lambda;
E_two = r01 + r12.*exp(1i*delta);
R_two = abs(E_two).^2;

% FFT
R_two_fft = abs(fft(R_two-mean(R_two))).^2;
freq = (0:N-1);  % FFT 样点序号

%% 情况2：多光束干涉（无限几何级数）
% 精确薄膜反射公式（Airy公式）
beta = 2*pi*n1*d*cos(theta1)./lambda; % 相位
r = (r01 + r12.*exp(2i*beta)) ./ (1 + r01.*r12.*exp(2i*beta));
R_multi = abs(r).^2;

% FFT
R_multi_fft = abs(fft(R_multi-mean(R_multi))).^2;

%% 绘图对比
figure;
subplot(2,2,1);
plot(nu,R_two,'r'); xlabel('波数 (cm^{-1})'); ylabel('R');
title('两束干涉反射率');

subplot(2,2,2);
plot(freq,R_two_fft,'b'); xlim([0 200]);
xlabel('FFT 样点'); ylabel('幅度');
title('两束干涉 FFT');

subplot(2,2,3);
plot(nu,R_multi,'r'); xlabel('波数 (cm^{-1})'); ylabel('R');
title('多光束干涉反射率');

subplot(2,2,4);
plot(freq,R_multi_fft,'b'); xlim([0 200]);
xlabel('FFT 样点'); ylabel('幅度');
title('多光束干涉 FFT');
