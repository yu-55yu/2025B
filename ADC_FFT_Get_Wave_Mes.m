function [VPP, Freq] = ADC_FFT_Get_Wave_Mes(Row, Fs, FFT_Output, correctNum, FFT_LEN)
    k = 2.667;
    DatePower1 = 0;
    DatePower2 = 0;
    
    % 在峰值周围±correctNum范围内计算能量重心
    for i = -correctNum:correctNum
        idx = Row + i;
        if idx >= 1 && idx <= length(FFT_Output)  % MATLAB索引从1开始，添加边界检查
            power = FFT_Output(idx) * FFT_Output(idx);
            DatePower1 = DatePower1 + idx * power;  % 用索引加权
            DatePower2 = DatePower2 + power;
        end
    end
    
    f = DatePower1 / DatePower2;  % 得到加权平均索引
    Freq = f * Fs / FFT_LEN;      % 校正后的频率
    VPP = 2.0 * sqrt(k * DatePower2);  % 幅值
end