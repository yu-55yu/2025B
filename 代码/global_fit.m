function [x_optimal, R_squared] = global_fit(x0, lb, ub, R_fit, waveLen_fit, angle, epsilon_inf, model_type)
    waveNum_fit = 10000 ./ waveLen_fit; 
    
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


