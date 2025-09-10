% -------------------------------------------------------------------
% �ű�����: ���ӻ�̼����(SiC)��ԲƬ�ڲ�ͬ������µķ������
%           - ��ʶ100%�����ʵ���������
%           - ��ǳ�������Ϊ0�����ݵ�
% ������Դ: ����1.xlsx (10��) �� ����2.xlsx (15��)
% -------------------------------------------------------------------

%% 1. ��ʼ������
clc;            % ��������д���
clear;          % ������������б���
close all;      % �ر�����ͼ�δ���

%% 2. �����ļ���Ϣ
file_sic_1.path = '����1.xlsx';
file_sic_1.angle = 10; % ����� (��)

file_sic_2.path = '����2.xlsx';
file_sic_2.angle = 15; % ����� (��)

%% 3. ��ȡ���ݲ��������
try
    % ��ȡ����1������
    data1 = readmatrix(file_sic_1.path);
    wavenumber1 = data1(:, 1);
    reflectance1 = data1(:, 2);
    % ���Ҹ���1�з�����Ϊ0�ĵ�
    zero_points1 = wavenumber1(reflectance1 == 0);

    % ��ȡ����2������
    data2 = readmatrix(file_sic_2.path);
    wavenumber2 = data2(:, 1);
    reflectance2 = data2(:, 2);
    % ���Ҹ���2�з�����Ϊ0�ĵ�
    zero_points2 = wavenumber2(reflectance2 == 0);

catch ME
    if (strcmp(ME.identifier, 'MATLAB:readmatrix:FileNotFound'))
        errorMessage = sprintf('����: δ�ҵ������ļ���\n��ȷ�� ''%s'' �� ''%s'' ��˽ű���ͬһĿ¼�¡�', file_sic_1.path, file_sic_2.path);
        uiwait(warndlg(errorMessage));
        return;
    else
        rethrow(ME);
    end
end

%% 4. ���ݿ��ӻ�
figure('Name', 'SiC Wafer Interference Spectrum', 'NumberTitle', 'off');
hold on;

% ����10���15������ǵ���������
plot(wavenumber1, reflectance1, 'LineWidth', 1.2, 'DisplayName', ['����� = ' num2str(file_sic_1.angle) '��']);
plot(wavenumber2, reflectance2, 'LineWidth', 1.2, 'DisplayName', ['����� = ' num2str(file_sic_2.angle) '��']);

% �� y=100 ��λ�û�һ����ɫ����
yline(100, '--r', 'LineWidth', 1.0, 'DisplayName', '100% ������');

% ������ʹ�� 'x' ��Ƿ�����Ϊ 0 �ĵ�
% ����Ƿ������㣬�����������б�ע
if ~isempty(zero_points1)
    scatter(zero_points1, zeros(size(zero_points1)), 40, 'kx', 'LineWidth', 1.5, 'DisplayName', '  0% ������');
end
% Ϊ����ͼ���ظ��������ı�ע�������DisplayName
if ~isempty(zero_points2)
    scatter(zero_points2, zeros(size(zero_points2)), 40, 'kx', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

hold off;

%% 5. ͼ�θ�ʽ��
xlabel('���� (cm^{-1})', 'FontSize', 12);
ylabel('������ (%)', 'FontSize', 12);
legend('show', 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

disp('���ݿ��ӻ���ɡ�');