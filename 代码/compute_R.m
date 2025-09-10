% function R = compute_R(thk, n1_complex, n2_complex, waveNum, theta0, opt)
% % compute_R 计算单层膜的反射率 (矢量化版本)
% %   thk: 膜厚 (μm)
% %   n1_complex: 外延层复折射率 (n+ik)
% %   n2_complex: 衬底复折射率 (n+ik)
% %   waveNum: 波数 (cm^-1)
% %   theta0: 入射角 (弧度)

%     n0 = 1.0;
%     n1 = n1_complex(:);
%     n2 = n2_complex(:);
%     waveNum = waveNum(:);
    
%     % Snell's Law
%     sin_theta1 = n0 * sin(theta0) ./ n1;
%     cos_theta1 = sqrt(1 - sin_theta1.^2);
%     sin_theta2 = n0 * sin(theta0) ./ n2;
%     cos_theta2 = sqrt(1 - sin_theta2.^2);
    
%     % Fresnel Coefficients
%     r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
%     r12_s = (n1.*cos_theta1 - n2.*cos_theta2) ./ (n1.*cos_theta1 + n2.*cos_theta2);
%     r01_p = (n1*cos(theta0) - n0*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
%     r12_p = (n2.*cos_theta1 - n1.*cos_theta2) ./ (n2.*cos_theta1 + n1.*cos_theta2);
    
%     % Phase shift
%     delta = 4 * pi * n1 .* thk .* cos_theta1 .* waveNum / 10000;
%     exp_term = exp(1i*delta);
    
%     if opt == 1
%         R = cal_R_mul(r01_s, r12_s, r01_p, r12_p, exp_term);
%     else 
%         R = cal_R_db(r01_s, r12_s, r01_p, r12_p, exp_term);
%     end
% end

function R = compute_R(thk, n1_complex, n2_complex, waveNum, theta0, opt)
% compute_R 计算单层膜的反射率
%   thk: 膜厚 (μm)
%   n1_complex: 外延层复折射率 (n+ik)
%   n2_complex: 衬底复折射率 (n+ik)
%   waveNum: 波数 (cm^-1)
%   theta0: 入射角 (弧度)
%   opt: 选项
%        0 - 强制使用双光束模型
%        1 - 强制使用多光束模型
%        2 - 自适应选择（推荐）

    n0 = 1.0;
    n1 = n1_complex(:);
    n2 = n2_complex(:);
    waveNum = waveNum(:);
    
    sin_theta1 = n0 * sin(theta0) ./ n1;
    cos_theta1 = sqrt(1 - sin_theta1.^2);
    sin_theta2 = n0 * sin(theta0) ./ n2;
    cos_theta2 = sqrt(1 - sin_theta2.^2);
    
    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    r12_s = (n1.*cos_theta1 - n2.*cos_theta2) ./ (n1.*cos_theta1 + n2.*cos_theta2);
    r01_p = (n1*cos(theta0) - n0*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    r12_p = (n2.*cos_theta1 - n1.*cos_theta2) ./ (n2.*cos_theta1 + n1.*cos_theta2);
    
    delta = 4 * pi * n1 .* thk .* cos_theta1 .* waveNum / 10000;
    exp_term = exp(1i*delta);
    
    if opt == 0
        % 强制双光束
        R = cal_R_db(r01_s, r12_s, r01_p, r12_p, exp_term);
    elseif opt == 1
        % 强制多光束
        R = cal_R_mul(r01_s, r12_s, r01_p, r12_p, exp_term);
    elseif opt == 2
        % 自适应选择模型
        epsilon_threshold = 0.01;
        [use_multi_beam, ~] = check(thk, n1_complex, n2_complex, waveNum, theta0, epsilon_threshold);
        R = zeros(size(waveNum));
        
        % 分别处理需要不同模型的波长点
        if any(use_multi_beam)
            idx = use_multi_beam;
            R(idx) = cal_R_mul(r01_s(idx), r12_s(idx), r01_p(idx), r12_p(idx), exp_term(idx));
        end
        if any(~use_multi_beam)
            idx = ~use_multi_beam;
            R(idx) = cal_R_db(r01_s(idx), r12_s(idx), r01_p(idx), r12_p(idx), exp_term(idx));
        end
        
    elseif opt >= 0.005 && opt <= 0.05
        % 使用opt作为自适应阈值
        epsilon_threshold = opt;
        [use_multi_beam, ~] = check(thk, n1_complex, n2_complex, waveNum, theta0, epsilon_threshold);
        R = zeros(size(waveNum));
        
        if any(use_multi_beam)
            idx = use_multi_beam;
            R(idx) = cal_R_mul(r01_s(idx), r12_s(idx), r01_p(idx), r12_p(idx), exp_term(idx));
        end
        if any(~use_multi_beam)
            idx = ~use_multi_beam;
            R(idx) = cal_R_db(r01_s(idx), r12_s(idx), r01_p(idx), r12_p(idx), exp_term(idx));
        end
    else
        % 默认使用多光束
        R = cal_R_mul(r01_s, r12_s, r01_p, r12_p, exp_term);
    end
end

