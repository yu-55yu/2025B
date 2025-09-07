function plot_fit_res(waveNum, R_exp, R_fit, angle, waveN_fit_min)
% PLOT_FIT_RES 绘制拟合结果与实测数据的对比图
%   waveNum: 波数 (cm^-1)
%   R_exp: 实验测得的反射率 (0-1)
%   R_fit: 模型拟合的反射率 (0-1)
%   angle: 入射角 (度)
%   waveN_fit_min: 拟合区域的起始波数

    color_blue = [0, 176, 232] / 255;
    color_red = [254, 100, 108] / 255;
    
    figure('Name', ['拟合结果分析 (入射角 ' num2str(angle) '°)'], 'Position', [150, 150, 1200, 600]);
    plot(waveNum, R_exp*100, '.', 'Color', color_blue, 'MarkerSize', 5, 'DisplayName', '实测数据');
    hold on;
    plot(waveNum, R_fit*100, '-', 'Color', color_red, 'LineWidth', 1.5, 'DisplayName', '拟合模型');
    
    ylim_vals = ylim;
    h = fill([waveN_fit_min, max(waveNum), max(waveNum), waveN_fit_min], ...
             [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
             'k', 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'DisplayName', '拟合区域');
    uistack(h, 'bottom');
    
    xlabel('波数 (cm^{-1})');
    ylabel('反射率 (%)');
    legend('Location', 'best');
    grid on;
    xlim([min(waveNum), max(waveNum)]);
end