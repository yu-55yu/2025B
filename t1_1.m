%% ̼�������Ӳ��Ȳⶨ������
clear; close all; clc;

%% 1. ��ȡ����
% ��ȡ����1�͸���2�����ݣ�10�Ⱥ�15������ǣ�
data1 = readmatrix('����1.xlsx'); % 10�������
data2 = readmatrix('����2.xlsx'); % 15�������

% ��ȡ�����ͷ�����
wavenumber1 = data1(:,1); % ���� (cm^-1)
reflectance1 = data1(:,2); % ������ (%)
wavenumber2 = data2(:,1);
reflectance2 = data2(:,2);

%% 2. ����Ԥ���� - ��ȡ��������1000������
idx1 = wavenumber1 >= 1000;
k1 = wavenumber1(idx1);
R1 = reflectance1(idx1);

idx2 = wavenumber2 >= 1000;
k2 = wavenumber2(idx2);
R2 = reflectance2(idx2);

%% 3. ʹ��������Ү������������Ӳ�������n1
% ������Ү������: n^2 = 1 + ��(Ai*��^2)/(��^2 - Bi)
% ���� �� = 10000/k (ת��Ϊ΢��)

% ѡ��������2000�����ݽ������
idx_fit = k1 > 2000;
k_fit = k1(idx_fit);
lambda_fit = 10000 ./ k_fit; % ת��Ϊ����(΢��)

% ̼���������ʵĳ�ʼ�ο�ֵ���������ף�
n_ref = 2.65; % SiC����������

% ����������Ү�����̵���Ϻ���
sellmeier = @(params, lambda) sqrt(1 + ...
    params(1)*lambda.^2./(lambda.^2 - params(4)) + ...
    params(2)*lambda.^2./(lambda.^2 - params(5)) + ...
    params(3)*lambda.^2./(lambda.^2 - params(6)));

% ��ʼ�������� [A1, A2, A3, B1, B2, B3]
initial_params = [5.0, 0.5, 0.1, 0.1, 0.5, 10];

% ʹ����С���˷����
% ��Ҫ�ӷ�������������ȡ��������Ϣ
% ���ڵ��㷴�䣬������ R = ((n1-n0)/(n1+n0))^2������n0=1(����)
% ��� n1 = (1+sqrt(R/100))/(1-sqrt(R/100))

% �����Ƶ����������ʣ�����ЧӦ��С��
R_high = R1(idx_fit);
n_estimated = (1 + sqrt(R_high/100)) ./ (1 - sqrt(R_high/100));

% ���������Ү������
options = optimset('Display','iter','TolFun',1e-8,'TolX',1e-8);
params_fit = lsqcurvefit(sellmeier, initial_params, lambda_fit, n_estimated, ...
    [0.1, 0.01, 0.001, 0.01, 0.1, 1], ...  % �½�
    [10, 5, 1, 1, 5, 50], ...               % �Ͻ�
    options);

% ��ȡ��ϲ���
A1 = params_fit(1); A2 = params_fit(2); A3 = params_fit(3);
B1 = params_fit(4); B2 = params_fit(5); B3 = params_fit(6);

fprintf('������Ү��������Ͻ��:\n');
fprintf('A1 = %.6f, A2 = %.6f, A3 = %.6f\n', A1, A2, A3);
fprintf('B1 = %.6f, B2 = %.6f, B3 = %.6f\n', B1, B2, B3);

% �������в�����Ӧ��������n1
lambda_all = 10000 ./ k1; % ���в���
n1_all = sellmeier(params_fit, lambda_all);

%% 4. FFT�����ҵ�������
% �Ӵ�����
win = hann(length(R1));
R1_windowed = R1 .* win;

% ����������
dk = mean(diff(k1)); % ƽ���������
Fs = 1 / dk;         % ����Ƶ��

% FFT����
N = length(R1_windowed);
Y = fft(R1_windowed);
P2 = abs(Y / N);
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2 * P1(2:end-1);

% Ƶ�����꣨��λ��cm��
f = Fs * (0:(N/2)) / N;

% �ҵ�����ֵ
threshold_f = 0.001; % Ƶ����ֵ��cm��
valid_idx = find(f > threshold_f);
[maxVal, peakIdx_temp] = max(P1(valid_idx));
peakIdx = valid_idx(peakIdx_temp);
f_peak = f(peakIdx);

% ��������У��������ȷ�ķ�ֵλ�ã�
correctNum = 5; % У����Χ
[f_corrected, ~] = energy_centroid_correction(peakIdx, f, P1, correctNum);

fprintf('\nԭʼ��Ƶ�� = %.4f cm\n', f_peak);
fprintf('��������У����Ƶ�� = %.4f cm\n', f_corrected);

%% 5. �ĵ�������n2�Ĺ���
% ʹ������Ƶ�ʴ�����������Ϊ�ĵ������ʵĳ�ʼ����
k_center = mean(k1);
lambda_center = 10000 / k_center;
n2_init = sellmeier(params_fit, lambda_center) * 0.95; % �ĵ�������ͨ����С

%% 6. ���Ӳ��ȵĳ�ʼ����
% ���ڴ�ֱ���䣬��̲� = 2*n1*d
% ����������2*n1*d = m*�ˣ�����m�Ǹ��漶��
% �ڲ�����2*n1*d*k = m*10000
% ���� ��k = 10000/(2*n1*d)
% ��� d = 10000/(2*n1*��k)

% ʹ��ƽ��������
n1_avg = mean(n1_all);
d_init = 1 / (2 * n1_avg * f_corrected); % ��ʼ��ȹ��ƣ�΢�ף�

fprintf('\n��ʼ��ȹ��� d = %.2f ��m\n', d_init * 10000);

%% 7. ��ȷ��ϣ���������Ǻ������ĸ���ģ��
theta1 = 10 * pi/180; % ����ǣ����ȣ�
theta2 = 15 * pi/180;

% ���淴����ģ�ͣ���������ǣ�
% R = |r1 + r2*exp(2i*beta)|^2 / |1 + r1*r2*exp(2i*beta)|^2
% ���� beta = 2*pi*n1*d*cos(theta_r)*k/10000
% theta_r������ǣ���Snell����ȷ��

reflectance_model = @(params, k, theta) calculate_reflectance(params, k, theta, params_fit);

% ������[d(΢��), n2, ��n(���Ӳ���ĵ������ʲ������)]
initial_params2 = [d_init*10000, n2_init, 0];

% ͬʱ��������Ƕȵ�����
k_combined = [k1; k2];
R_combined = [R1; R2];
theta_combined = [ones(size(k1))*theta1; ones(size(k2))*theta2];

% ��������С�������
options2 = optimoptions('lsqcurvefit','Display','iter','MaxFunctionEvaluations',5000);
params_final = lsqcurvefit(@(p,k) reflectance_model(p, k_combined, theta_combined), ...
    initial_params2, k_combined, R_combined, ...
    [d_init*10000*0.5, n2_init*0.9, -0.1], ...  % �½�
    [d_init*10000*1.5, n2_init*1.1, 0.1], ...   % �Ͻ�
    options2);

d_final = params_final(1);      % ��ȣ�΢�ף�
n2_final = params_final(2);     % �ĵ�������
delta_n = params_final(3);      % ����������

fprintf('\n������Ͻ��:\n');
fprintf('���Ӳ��� d = %.3f ��m\n', d_final);
fprintf('�ĵ������� n2 = %.4f\n', n2_final);
fprintf('���������� ��n = %.4f\n', delta_n);

%% 8. ������ӻ�
figure('Position', [100, 100, 1200, 800]);

% ��ͼ1��ʵ������϶Աȣ�10�ȣ�
subplot(2,2,1);
R_fit1 = reflectance_model(params_final, k1, theta1*ones(size(k1)));
plot(k1, R1, 'b-', 'LineWidth', 1, 'DisplayName', 'ʵ������');
hold on;
plot(k1, R_fit1, 'r--', 'LineWidth', 1, 'DisplayName', '��Ͻ��');
xlabel('���� (cm^{-1})');
ylabel('������ (%)');
title('10�������');
legend('Location', 'best');
grid on;

% ��ͼ2��ʵ������϶Աȣ�15�ȣ�
subplot(2,2,2);
R_fit2 = reflectance_model(params_final, k2, theta2*ones(size(k2)));
plot(k2, R2, 'b-', 'LineWidth', 1, 'DisplayName', 'ʵ������');
hold on;
plot(k2, R_fit2, 'r--', 'LineWidth', 1, 'DisplayName', '��Ͻ��');
xlabel('���� (cm^{-1})');
ylabel('������ (%)');
title('15�������');
legend('Location', 'best');
grid on;

% ��ͼ3��FFTƵ��
subplot(2,2,3);
plot(f, P1, 'b-', 'LineWidth', 1);
hold on;
plot(f_peak, P1(peakIdx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
plot(f_corrected, P1(peakIdx), 'g^', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('�ռ�Ƶ�� (cm)');
ylabel('����');
title('FFTƵ�׷���');
legend('FFTƵ��', 'ԭʼ��ֵ', 'У�����ֵ', 'Location', 'best');
xlim([0, 0.01]);
grid on;

% ��ͼ4���������沨���仯
subplot(2,2,4);
plot(lambda_all, n1_all, 'b-', 'LineWidth', 2);
hold on;
yline(n2_final, 'r--', 'LineWidth', 1.5);
xlabel('���� (��m)');
ylabel('������');
title('�����ʷֲ�');
legend('���Ӳ� n_1', '�ĵ� n_2', 'Location', 'best');
grid on;

%% 9. ������
% ������ϲв�
residual1 = R1 - R_fit1;
residual2 = R2 - R_fit2;
RMSE1 = sqrt(mean(residual1.^2));
RMSE2 = sqrt(mean(residual2.^2));

fprintf('\n���������:\n');
fprintf('10������� RMSE = %.4f%%\n', RMSE1);
fprintf('15������� RMSE = %.4f%%\n', RMSE2);

%% ��������
function R = calculate_reflectance(params, k, theta, sellmeier_params)
    % ���㷴����
    d = params(1);      % ��ȣ�΢�ף�
    n2 = params(2);     % �ĵ�������
    delta_n = params(3); % ����������
    
    % �������Ӳ�������
    lambda = 10000 ./ k;
    n1 = sqrt(1 + ...
        sellmeier_params(1)*lambda.^2./(lambda.^2 - sellmeier_params(4)) + ...
        sellmeier_params(2)*lambda.^2./(lambda.^2 - sellmeier_params(5)) + ...
        sellmeier_params(3)*lambda.^2./(lambda.^2 - sellmeier_params(6))) + delta_n;
    
    % ����������
    n0 = 1;
    
    % Snell���ɼ��������
    sin_theta_r = sin(theta) .* n0 ./ n1;
    cos_theta_r = sqrt(1 - sin_theta_r.^2);
    
    % ����������ϵ����sƫ��
    r01 = (n0*cos(theta) - n1.*cos_theta_r) ./ (n0*cos(theta) + n1.*cos_theta_r);
    r12 = (n1.*cos_theta_r - n2.*cos_theta_r) ./ (n1.*cos_theta_r + n2.*cos_theta_r);
    
    % ��λ
    beta = 2*pi*n1.*d.*cos_theta_r.*k/10000;
    
    % �ܷ�����
    numerator = abs(r01 + r12.*exp(2i*beta)).^2;
    denominator = abs(1 + r01.*r12.*exp(2i*beta)).^2;
    R = 100 * numerator ./ denominator;
end

function [f_corrected, P_corrected] = energy_centroid_correction(peakIdx, f, P, correctNum)
    % ��������У��
    idx_range = max(1, peakIdx-correctNum):min(length(f), peakIdx+correctNum);
    f_range = f(idx_range);
    P_range = P(idx_range);
    
    % ������������
    f_corrected = sum(f_range .* P_range) / sum(P_range);
    P_corrected = interp1(f, P, f_corrected);
end