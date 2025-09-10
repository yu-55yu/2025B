% -------------------------------------------------------------------
% 脚本功能: 可视化碳化硅(SiC)晶圆片在不同入射角下的反射光谱
%           - 标识100%反射率的理论上限
%           - 标记出反射率为0的数据点
% 数据来源: 附件1.xlsx (10°) 和 附件2.xlsx (15°)
% -------------------------------------------------------------------

%% 1. 初始化环境
clc;            % 清空命令行窗口
clear;          % 清除工作区所有变量
close all;      % 关闭所有图形窗口

%% 2. 定义文件信息
file_sic_1.path = '附件1.xlsx';
file_sic_1.angle = 10; % 入射角 (度)

file_sic_2.path = '附件2.xlsx';
file_sic_2.angle = 15; % 入射角 (度)

%% 3. 读取数据并查找零点
try
    % 读取附件1的数据
    data1 = readmatrix(file_sic_1.path);
    wavenumber1 = data1(:, 1);
    reflectance1 = data1(:, 2);
    % 查找附件1中反射率为0的点
    zero_points1 = wavenumber1(reflectance1 == 0);

    % 读取附件2的数据
    data2 = readmatrix(file_sic_2.path);
    wavenumber2 = data2(:, 1);
    reflectance2 = data2(:, 2);
    % 查找附件2中反射率为0的点
    zero_points2 = wavenumber2(reflectance2 == 0);

catch ME
    if (strcmp(ME.identifier, 'MATLAB:readmatrix:FileNotFound'))
        errorMessage = sprintf('错误: 未找到数据文件。\n请确保 ''%s'' 和 ''%s'' 与此脚本在同一目录下。', file_sic_1.path, file_sic_2.path);
        uiwait(warndlg(errorMessage));
        return;
    else
        rethrow(ME);
    end
end

%% 4. 数据可视化
figure('Name', 'SiC Wafer Interference Spectrum', 'NumberTitle', 'off');
hold on;

% 绘制10°和15°入射角的数据曲线
plot(wavenumber1, reflectance1, 'LineWidth', 1.2, 'DisplayName', ['入射角 = ' num2str(file_sic_1.angle) '°']);
plot(wavenumber2, reflectance2, 'LineWidth', 1.2, 'DisplayName', ['入射角 = ' num2str(file_sic_2.angle) '°']);

% 在 y=100 的位置画一条红色虚线
yline(100, '--r', 'LineWidth', 1.0, 'DisplayName', '100% 反射率');

% 新增：使用 'x' 标记反射率为 0 的点
% 检查是否存在零点，如果存在则进行标注
if ~isempty(zero_points1)
    scatter(zero_points1, zeros(size(zero_points1)), 40, 'kx', 'LineWidth', 1.5, 'DisplayName', '  0% 反射率');
end
% 为避免图例重复，后续的标注不再添加DisplayName
if ~isempty(zero_points2)
    scatter(zero_points2, zeros(size(zero_points2)), 40, 'kx', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

hold off;

%% 5. 图形格式化
xlabel('波数 (cm^{-1})', 'FontSize', 12);
ylabel('反射率 (%)', 'FontSize', 12);
legend('show', 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

disp('数据可视化完成。');