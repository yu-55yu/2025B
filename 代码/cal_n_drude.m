function n2_complex = cal_n_drude(waveNum, drudeParam, epsilon_inf)
    nu_p = drudeParam(1);
    Gamma = drudeParam(2);
    nu = waveNum(:);
    epsilon_complex = epsilon_inf - (nu_p^2) ./ (nu.^2 + 1i * Gamma * nu);
    n2_complex = sqrt(epsilon_complex);
end