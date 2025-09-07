function n1_complex = cal_n_sellmeier(waveLen, selParam, k1Param)
    B1=selParam(1); B2=selParam(2); B3=selParam(3);
    C1=selParam(4); C2=selParam(5); C3=selParam(6);
    lambda_sq = waveLen.^2;
    eps = 1e-10;
    term1 = B1*lambda_sq./(lambda_sq-C1+eps);
    term2 = B2*lambda_sq./(lambda_sq-C2+eps);
    term3 = B3*lambda_sq./(lambda_sq-C3+eps);
    n_squared = 1 + term1 + term2 + term3;
    n_squared(n_squared < 1) = 1;
    n_real = real(sqrt(n_squared));
    k1_A = k1Param(1); k1_B = k1Param(2);
    k_imag = k1_A * waveLen.^k1_B;
    n1_complex = n_real + 1i * k_imag;
end