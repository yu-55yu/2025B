% function R = compute_R(thk, n1_complex, n2_complex, waveNum, theta0, opt)
% % compute_R ���㵥��Ĥ�ķ����� (ʸ�����汾)
% %   thk: Ĥ�� (��m)
% %   n1_complex: ���Ӳ㸴������ (n+ik)
% %   n2_complex: �ĵ׸������� (n+ik)
% %   waveNum: ���� (cm^-1)
% %   theta0: ����� (����)

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
% compute_R ���㵥��Ĥ�ķ�����
%   thk: Ĥ�� (��m)
%   n1_complex: ���Ӳ㸴������ (n+ik)
%   n2_complex: �ĵ׸������� (n+ik)
%   waveNum: ���� (cm^-1)
%   theta0: ����� (����)
%   opt: ѡ��
%        0 - ǿ��ʹ��˫����ģ��
%        1 - ǿ��ʹ�ö����ģ��
%        2 - ����Ӧѡ���Ƽ���

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
        % ǿ��˫����
        R = cal_R_db(r01_s, r12_s, r01_p, r12_p, exp_term);
    elseif opt == 1
        % ǿ�ƶ����
        R = cal_R_mul(r01_s, r12_s, r01_p, r12_p, exp_term);
    elseif opt == 2
        % ����Ӧѡ��ģ��
        epsilon_threshold = 0.01;
        [use_multi_beam, ~] = check(thk, n1_complex, n2_complex, waveNum, theta0, epsilon_threshold);
        R = zeros(size(waveNum));
        
        % �ֱ�����Ҫ��ͬģ�͵Ĳ�����
        if any(use_multi_beam)
            idx = use_multi_beam;
            R(idx) = cal_R_mul(r01_s(idx), r12_s(idx), r01_p(idx), r12_p(idx), exp_term(idx));
        end
        if any(~use_multi_beam)
            idx = ~use_multi_beam;
            R(idx) = cal_R_db(r01_s(idx), r12_s(idx), r01_p(idx), r12_p(idx), exp_term(idx));
        end
        
    elseif opt >= 0.005 && opt <= 0.05
        % ʹ��opt��Ϊ����Ӧ��ֵ
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
        % Ĭ��ʹ�ö����
        R = cal_R_mul(r01_s, r12_s, r01_p, r12_p, exp_term);
    end
end

