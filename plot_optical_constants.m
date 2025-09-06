function plot_optical_constants(waveLen, waveNum, n_complex, layer_name)
% PLOT_OPTICAL_CONSTANTS 绘制复折射率色散曲线
%   waveLen: 波长 (μm)
%   waveNum: 波数 (cm^-1)
%   n_complex: 复折射率 (n+ik)
%   layer_name: 图例中显示的层名称 (e.g., '外延层')

    figure('Name', ['拟合得到的 ' layer_name ' 复折射率色散曲线'], 'Position', [100, 100, 1200, 450]);
    
    n_real = real(n_complex);
    k_imag = imag(n_complex);
    
    color_n = [0, 114, 189] / 255; % 蓝色
    color_k = [217, 83, 25] / 255;  % 红色
    
    % 图1: vs 波数
    ax1 = nexttile;
    hold(ax1, 'on');
    yyaxis(ax1, 'left');
    plot(ax1, waveNum, n_real, '-', 'Color', color_n, 'LineWidth', 1.5);
    ylabel(ax1, '折射率实部 n');
    ax1.YColor = color_n;
    
    yyaxis(ax1, 'right');
    plot(ax1, waveNum, k_imag, '--', 'Color', color_k, 'LineWidth', 1.5);
    ylabel(ax1, '消光系数 k');
    ax1.YColor = color_k;
    set(ax1, 'YScale', 'log'); % k值通常用对数坐标
    
    hold(ax1, 'off');
    grid(ax1, 'on');
    xlabel(ax1, '波数 (cm^{-1})');
    xlim(ax1, [min(waveNum), max(waveNum)]);

    % 图2: vs 波长
    ax2 = nexttile;
    hold(ax2, 'on');
    yyaxis(ax2, 'left');
    plot(ax2, waveLen, n_real, '-', 'Color', color_n, 'LineWidth', 1.5);
    ylabel(ax2, '折射率实部 n');
    ax2.YColor = color_n;
    
    yyaxis(ax2, 'right');
    plot(ax2, waveLen, k_imag, '--', 'Color', color_k, 'LineWidth', 1.5);
    ylabel(ax2, '消光系数 k');
    ax2.YColor = color_k;
    set(ax2, 'YScale', 'log');
    
    hold(ax2, 'off');
    grid(ax2, 'on');
    xlabel(ax2, '波长 (μm)');
    xlim(ax2, [min(waveLen), max(waveLen)]);
end