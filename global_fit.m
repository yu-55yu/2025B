function [x_optimal, R_squared] = global_fit(x0, lb, ub, waveNum_fit, R_fit, waveLen_fit, angle, epsilon_inf, model_type)
    theta0_rad = angle * pi / 180;
    if strcmp(model_type, 'cauchy')
        model_func = @(x, k) model_R_cauchy(x, k, waveLen_fit, theta0_rad, epsilon_inf);
    elseif strcmp(model_type, 'sellmeier')
        model_func = @(x, k) model_R_sellmeier(x, k, waveLen_fit, theta0_rad, epsilon_inf);
    elseif strcmp(model_type, 'real_substrate')
        model_func = @(x, k) model_R_real_substrate(x, k, waveLen_fit, theta0_rad);
    end
    options = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxIterations', 1000, 'FunctionTolerance', 1e-9, 'StepTolerance', 1e-10, 'Algorithm', 'trust-region-reflective', 'UseParallel', true);
    [x_optimal, resnorm] = lsqcurvefit(model_func, x0, waveNum_fit, R_fit, lb, ub, options);
    SS_tot = sum((R_fit - mean(R_fit)).^2);
    R_squared = 1 - resnorm / SS_tot;
end

function R = model_R_sellmeier(params, ~, waveLen, theta0, epsilon_inf)
    thk = params(1);
    selParam = params(2:7);
    k1Param = params(8:9);
    drudeParam = params(10:11);
    waveNum_model = 10000 ./ waveLen;
    n1_complex = cal_n_sellmeier(waveLen, selParam, k1Param);
    n2_complex = cal_n_drude(waveNum_model, drudeParam, epsilon_inf);
    R = compute_R(thk, n1_complex, n2_complex, waveNum_model, theta0,1);
end


function R = model_R_real_substrate(params, ~, waveLen, theta0)
    thk = params(1);
    n2 = params(2);
    selParam = params(3:8);
    k1Param = params(9:10);
    waveNum_model = 10000 ./ waveLen;
    n1_complex = cal_n_sellmeier(waveLen, selParam, k1Param);
    R = compute_R(thk, n1_complex, n2, waveNum_model, theta0, 1);
end



function R = model_R_cauchy(params, ~, waveLen, theta0, epsilon_inf)
    thk = params(1);
    cauchyParam = params(2:3);
    k1Param = params(4:5);
    drudeParam = params(6:7);
    waveNum_model = 10000 ./ waveLen;
    n1_complex = cal_n_cauchy(waveLen, cauchyParam, k1Param);
    n2_complex = cal_n_drude(waveNum_model, drudeParam, epsilon_inf);
    R = compute_R(thk, n1_complex, n2_complex, waveNum_model, theta0,1);
end



function n1_complex = cal_n_cauchy(waveLen, cauchyParam, k1Param)
    A = cauchyParam(1); B = cauchyParam(2);
    n_real = A + B ./ (waveLen.^2);
    k1_A = k1Param(1); k1_B = k1Param(2);
    k_imag = k1_A * waveLen.^k1_B;
    n1_complex = n_real + 1i * k_imag;
end

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

function n2_complex = cal_n_drude(waveNum, drudeParam, epsilon_inf)
    nu_p = drudeParam(1);
    Gamma = drudeParam(2);
    nu = waveNum(:);
    epsilon_complex = epsilon_inf - (nu_p^2) ./ (nu.^2 + 1i * Gamma * nu);
    n2_complex = sqrt(epsilon_complex);
end