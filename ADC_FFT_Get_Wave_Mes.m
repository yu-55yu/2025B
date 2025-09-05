function [VPP, Freq] = ADC_FFT_Get_Wave_Mes(Row, Fs, FFT_Output, correctNum, FFT_LEN)
    k = 2.667;
    DatePower1 = 0;
    DatePower2 = 0;
    
    for i = -correctNum:correctNum
        idx = Row + i;
        if idx >= 1 && idx <= length(FFT_Output) 
            power = FFT_Output(idx) * FFT_Output(idx);
            DatePower1 = DatePower1 + idx * power;
            DatePower2 = DatePower2 + power;
        end
    end
    
    f = DatePower1 / DatePower2;
    Freq = f * Fs / FFT_LEN; 
    VPP = 2.0 * sqrt(k * DatePower2); 
end