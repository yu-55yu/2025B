function n1_complex = cal_n_cauchy(waveLen, cauchyParam, k1Param)
    A = cauchyParam(1); B = cauchyParam(2);
    n_real = A + B ./ (waveLen.^2);
    k1_A = k1Param(1); k1_B = k1Param(2);
    k_imag = k1_A * waveLen.^k1_B;
    n1_complex = n_real + 1i * k_imag;
end