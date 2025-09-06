clear; % �������������
clc; % ��������д���
close all; % �ر�����ͼ�δ���

filename = '����2.xlsx';
 dataTable = readtable(filename);

x = dataTable{:, 1}; % ��ȡ��һ��: ���� (cm^-1)
y = dataTable{:, 2}; % ��ȡ�ڶ���: ������ (%)

figure;
plot(x, y, 'b-', 'LineWidth', 1.5);
title('ԭʼ����: ������ vs ����');
xlabel('���� (cm^{-1})');
ylabel('������ (%)');
grid on;
set(gca, 'FontSize', 12);

% idx = (x>=400)&(x<730);
idx = x >= 1000; % �ҵ��±� �� 1000 cm^-1 ��λ��
x_cut = x(idx); % ��ȡ��Ĳ���
y_cut = y(idx); % ��ȡ��ķ�����

win = hann(length(y_cut));
y_win = y_cut .* win; 

dx = mean(diff(x_cut)); % ƽ���������(cm^-1)
Fs = 1 / dx; % "����Ƶ��"
N = length(y_win); % �źų���
Y = fft(y_win); % �Ӵ���� FFT
P2 = abs(Y / N); % ˫�����
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2 * P1(2:end-1);
% �ռ������꣨��λ cm��
f = Fs * (0:(N/2)) / N;


% �ҵ�����ֵ
threshold_f = 0.001;           % Ƶ����ֵ (cm)
valid_idx = find(f > threshold_f); % �ҵ����� f > 0.01 cm ������

% --- �� f > 0.01 cm ��Χ�������ֵ ---
[maxVal, peakIdx] = max(P1(valid_idx)); % ע�⣺peakIdx �� valid_idx �ľֲ�����
peakIdx = valid_idx(peakIdx); % ת��Ϊȫ������
f_peak = f(peakIdx); % ԭʼ��ֵƵ��

% ��������У��
correctNum = 5; % У����Χ (�ɵ�)
[VPP_corrected, f_corrected] = ADC_FFT_Get_Wave_Mes(peakIdx, Fs, P1, correctNum, N);

% ��ʾ���
disp(['ԭʼ����Ƶ�� = ', num2str(f_peak), ' cm']);
disp(['��������У��Ƶ�� = ', num2str(f_corrected), ' cm']);
% disp(['У������ֵ = ', num2str(VPP_corrected)]);


figure;
plot(f, P1, 'r-', 'LineWidth', 1.5);
hold on;
xlabel('�ռ�λ�� (cm)');
ylabel('���');

grid on;
set(gca, 'FontSize', 12);
xlim([0, max(f)]);
