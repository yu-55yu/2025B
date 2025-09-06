function plot_optical_constants(waveLen, waveNum, n_complex, layer_name)
% PLOT_OPTICAL_CONSTANTS ���Ƹ�������ɫɢ����
%   waveLen: ���� (��m)
%   waveNum: ���� (cm^-1)
%   n_complex: �������� (n+ik)
%   layer_name: ͼ������ʾ�Ĳ����� (e.g., '���Ӳ�')

    figure('Name', ['��ϵõ��� ' layer_name ' ��������ɫɢ����'], 'Position', [100, 100, 1200, 450]);
    
    n_real = real(n_complex);
    k_imag = imag(n_complex);
    
    color_n = [0, 114, 189] / 255; % ��ɫ
    color_k = [217, 83, 25] / 255;  % ��ɫ
    
    % ͼ1: vs ����
    ax1 = nexttile;
    hold(ax1, 'on');
    yyaxis(ax1, 'left');
    plot(ax1, waveNum, n_real, '-', 'Color', color_n, 'LineWidth', 1.5);
    ylabel(ax1, '������ʵ�� n');
    ax1.YColor = color_n;
    
    yyaxis(ax1, 'right');
    plot(ax1, waveNum, k_imag, '--', 'Color', color_k, 'LineWidth', 1.5);
    ylabel(ax1, '����ϵ�� k');
    ax1.YColor = color_k;
    set(ax1, 'YScale', 'log'); % kֵͨ���ö�������
    
    hold(ax1, 'off');
    grid(ax1, 'on');
    xlabel(ax1, '���� (cm^{-1})');
    xlim(ax1, [min(waveNum), max(waveNum)]);

    % ͼ2: vs ����
    ax2 = nexttile;
    hold(ax2, 'on');
    yyaxis(ax2, 'left');
    plot(ax2, waveLen, n_real, '-', 'Color', color_n, 'LineWidth', 1.5);
    ylabel(ax2, '������ʵ�� n');
    ax2.YColor = color_n;
    
    yyaxis(ax2, 'right');
    plot(ax2, waveLen, k_imag, '--', 'Color', color_k, 'LineWidth', 1.5);
    ylabel(ax2, '����ϵ�� k');
    ax2.YColor = color_k;
    set(ax2, 'YScale', 'log');
    
    hold(ax2, 'off');
    grid(ax2, 'on');
    xlabel(ax2, '���� (��m)');
    xlim(ax2, [min(waveLen), max(waveLen)]);
end