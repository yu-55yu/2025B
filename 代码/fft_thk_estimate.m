function thk = fft_thk_estimate(waveNum, R, n_avg, waveN_fit_min, theta0_deg)
% FFT_THK_ESTIMATE 通过傅里叶变换估算薄膜厚度
%   waveNum: 波数 (cm^-1)
%   R: 反射率
%   n_avg: 外延层在该区域的平均折射率
%   waveN_fit_min: 用于FFT的最小波数
%   theta0_deg: 入射角 (度)

    if nargin < 5, theta0_deg = 0; end
    theta0_rad = theta0_deg * pi / 180;
    cos_theta1 = real(sqrt(1 - (sin(theta0_rad) / n_avg)^2));
    
    filter = waveNum > waveN_fit_min;
    waveNum_fft = waveNum(filter);
    R_fft = R(filter);
    
    % 移除直流分量
    R_ac = R_fft - mean(R_fft);
    
    % 插值到均匀间隔的波数点
    N = 2^nextpow2(8*length(waveNum_fft));
    k_uniform = linspace(min(waveNum_fft), max(waveNum_fft), N);
    r_uniform = interp1(waveNum_fft, R_ac, k_uniform, 'pchip', 'extrap');
    
    % 加窗并执行FFT
    window = hann(N);
    r_win = r_uniform(:) .* window(:);
    fft_result = fft(r_win);
    fft_power = abs(fft_result(1:N/2)).^2;
    
    % 计算厚度轴
    dk = mean(diff(k_uniform));
    thk_axis = (0:N/2-1) * 10000 / (2 * n_avg * cos_theta1 * N * dk);
    
    % 寻找主峰
    search_range = find(thk_axis > 5 & thk_axis < 200); % 在合理范围内搜索
    if isempty(search_range)
        thk = 20; % 如果找不到，返回一个默认值
        return;
    end
    
    [~, max_idx_in_range] = max(fft_power(search_range));
    max_idx_rough = search_range(max_idx_in_range);
    
    % 重心法修正峰位
    correctNum = 3;
    DatePower1 = 0;
    DatePower2 = 0;
    for i = -correctNum:correctNum
        idx = max_idx_rough + i;
        if idx >= 1 && idx <= length(fft_power)
            power = fft_power(idx);
            DatePower1 = DatePower1 + idx * power;
            DatePower2 = DatePower2 + power;
        end
    end
    
    if DatePower2 > 0
        f_corrected = DatePower1 / DatePower2;
    else
        f_corrected = max_idx_rough;
    end
    
    thk = (f_corrected - 1) * 10000 / (2 * n_avg * cos_theta1 * N * dk);
end