function plot_fit_res(waveNum, R_exp, R_fit, angle, waveN_fit_min)
% PLOT_FIT_RES ������Ͻ����ʵ�����ݵĶԱ�ͼ
%   waveNum: ���� (cm^-1)
%   R_exp: ʵ���õķ����� (0-1)
%   R_fit: ģ����ϵķ����� (0-1)
%   angle: ����� (��)
%   waveN_fit_min: ����������ʼ����

    color_blue = [0, 176, 232] / 255;
    color_red = [254, 100, 108] / 255;
    
    figure('Name', ['��Ͻ������ (����� ' num2str(angle) '��)'], 'Position', [150, 150, 1200, 600]);
    plot(waveNum, R_exp*100, '.', 'Color', color_blue, 'MarkerSize', 5, 'DisplayName', 'ʵ������');
    hold on;
    plot(waveNum, R_fit*100, '-', 'Color', color_red, 'LineWidth', 1.5, 'DisplayName', '���ģ��');
    
    ylim_vals = ylim;
    h = fill([waveN_fit_min, max(waveNum), max(waveNum), waveN_fit_min], ...
             [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
             'k', 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'DisplayName', '�������');
    uistack(h, 'bottom');
    
    xlabel('���� (cm^{-1})');
    ylabel('������ (%)');
    legend('Location', 'best');
    grid on;
    xlim([min(waveNum), max(waveNum)]);
end