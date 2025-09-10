function [use_multi_beam, q_values] = check(thk, n1_complex, n2_complex, waveNum, theta0, epsilon_threshold)
% CHECK_MULTI_BEAM_INTERFERENCE 判断是否发生多光束干涉
%   根据判据 |q| = |r10*r12| * exp(-sigma*2*d) >= epsilon 来判断
%
% 输入参数:
%   thk: 膜厚 (μm)
%   n1_complex: 外延层复折射率 (n+ik)
%   n2_complex: 衬底复折射率 (n+ik)
%   waveNum: 波数 (cm^-1)
%   theta0: 入射角 (弧度)
%   epsilon_threshold: 阈值，默认0.03，通常取[0.01, 0.05]
%
% 输出:
%   use_multi_beam: 逻辑数组，true表示该波长需要考虑多光束干涉
%   q_values: 各波长对应的q值，用于调试和分析

    % 设置默认阈值
    if nargin < 6
        epsilon_threshold = 0.03;  % 默认取中间值
    end
    
    % 确保输入为列向量
    n0 = 1.0;
    n1 = n1_complex(:);
    n2 = n2_complex(:);
    waveNum = waveNum(:);
    
    % 计算波长 (μm)
    waveLen = 10000 ./ waveNum;
    
    % 提取复折射率的实部和虚部
    n1_real = real(n1);
    k1_imag = imag(n1);
    n2_real = real(n2);
    k2_imag = imag(n2);
    
    % 应用Snell定律计算折射角
    sin_theta1 = n0 * sin(theta0) ./ n1;
    cos_theta1 = sqrt(1 - sin_theta1.^2);
    sin_theta2 = n0 * sin(theta0) ./ n2;
    cos_theta2 = sqrt(1 - sin_theta2.^2);
    
    % 计算Fresnel反射系数 (s偏振和p偏振)
    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    r12_s = (n1.*cos_theta1 - n2.*cos_theta2) ./ (n1.*cos_theta1 + n2.*cos_theta2);
    r01_p = (n1*cos(theta0) - n0*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    r12_p = (n2.*cos_theta1 - n1.*cos_theta2) ./ (n2.*cos_theta1 + n1.*cos_theta2);
    
    % 计算平均反射系数的模
    r01_avg = sqrt((abs(r01_s).^2 + abs(r01_p).^2) / 2);
    r12_avg = sqrt((abs(r12_s).^2 + abs(r12_p).^2) / 2);
    
    % 计算吸收系数 α = 4πk/λ (单位: μm^-1)
    alpha = 4 * pi * k1_imag ./ waveLen;
    
    % 计算σ = 2α (考虑斜入射时的实际光程)
    % 在薄膜中的实际传播距离需要考虑折射角
    actual_path_factor = real(cos_theta1);  % 实际光程修正因子
    sigma = 2 * alpha ./ actual_path_factor;
    
    % 计算衰减因子 exp(-σ*2d)
    % 2d表示光在薄膜中往返一次的距离
    attenuation_factor = exp(-sigma * 2 * thk);
    
    % 计算判据值 |q| = |r10*r12| * exp(-σ*2d)
    q_values = r01_avg .* r12_avg .* attenuation_factor;



figure; 
plot(waveNum, q_values, 'LineWidth', 1.2);
hold on;


e_r = [0.6350, 0.0780, 0.1840]; % 这是 MATLAB 的  色

yline(epsilon_threshold, 'LineStyle', '--', 'Color',e_r, 'LineWidth', 1.2);

hold off; % 操作完成，关闭 hold

xlabel('波数 (cm^{-1})');
ylabel('q 值');
grid on; 
xlim([min(waveNum), max(waveNum)]);

legend('q 值', sprintf('阈值 ε=%.2f', epsilon_threshold));
    
    % 判断是否需要考虑多光束干涉
    use_multi_beam = q_values >= epsilon_threshold;
    
    % 输出统计信息（可选）
    if nargout ~= 0
        fprintf('\n=== 多光束干涉判断结果 ===\n');
        fprintf('阈值 ε = %.3f\n', epsilon_threshold);
        fprintf('q值范围: [%.4f, %.4f]\n', min(q_values), max(q_values));
        fprintf('需要多光束模型的波长点: %d/%d (%.1f%%)\n', ...
            sum(use_multi_beam), length(use_multi_beam), ...
            100*sum(use_multi_beam)/length(use_multi_beam));
        
        % 找出临界波长区域
        transition_idx = find(diff(use_multi_beam) ~= 0);
        if ~isempty(transition_idx)
            fprintf('模型转换波长点 (cm^-1): ');
            fprintf('%.1f ', waveNum(transition_idx));
            fprintf('\n');
        end
    end
end

function R = cal_R_mul(r01_s, r12_s, r01_p, r12_p, exp_term)
% CALCULATE_R_MULTI_BEAM 多光束干涉模型（完整的Airy公式）
    % s偏振
    numerator_s = abs(r01_s + r12_s .* exp_term).^2;
    denominator_s = abs(1 + r01_s .* r12_s .* exp_term).^2;
    R_s = numerator_s ./ denominator_s;
    
    % p偏振
    numerator_p = abs(r01_p + r12_p .* exp_term).^2;
    denominator_p = abs(1 + r01_p .* r12_p .* exp_term).^2;
    R_p = numerator_p ./ denominator_p;
    
    % 非偏振光的反射率（s和p的平均）
    R = real((R_s + R_p) / 2);
end

function R = compute_R_adaptive(thk, n1_complex, n2_complex, waveNum, theta0, epsilon_threshold)
% COMPUTE_R_ADAPTIVE 自适应选择计算模型的反射率函数
%   根据物理条件自动选择双光束或多光束干涉模型
%
% 输入参数：
%   thk: 膜厚 (μm)
%   n1_complex: 外延层复折射率
%   n2_complex: 衬底复折射率
%   waveNum: 波数 (cm^-1)
%   theta0: 入射角 (弧度)
%   epsilon_threshold: 多光束判据阈值 (默认0.03)

    if nargin < 6
        epsilon_threshold = 0.03;
    end
    
    % 判断是否需要多光束模型
    [use_multi_beam, q_values] = check_multi_beam_interference(thk, n1_complex, n2_complex, waveNum, theta0, epsilon_threshold);
    
    % 准备计算所需的参数
    n0 = 1.0;
    n1 = n1_complex(:);
    n2 = n2_complex(:);
    waveNum = waveNum(:);
    
    % Snell's Law
    sin_theta1 = n0 * sin(theta0) ./ n1;
    cos_theta1 = sqrt(1 - sin_theta1.^2);
    sin_theta2 = n0 * sin(theta0) ./ n2;
    cos_theta2 = sqrt(1 - sin_theta2.^2);
    
    % Fresnel Coefficients
    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    r12_s = (n1.*cos_theta1 - n2.*cos_theta2) ./ (n1.*cos_theta1 + n2.*cos_theta2);
    r01_p = (n1*cos(theta0) - n0*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    r12_p = (n2.*cos_theta1 - n1.*cos_theta2) ./ (n2.*cos_theta1 + n1.*cos_theta2);
    
    % Phase shift
    delta = 4 * pi * n1 .* thk .* cos_theta1 .* waveNum / 10000;
    exp_term = exp(1i*delta);
    
    % 初始化结果数组
    R = zeros(size(waveNum));
    
    % 对需要多光束模型的波长点使用多光束计算
    multi_beam_idx = use_multi_beam;
    if any(multi_beam_idx)
        R(multi_beam_idx) = cal_R_mul(...
            r01_s(multi_beam_idx), r12_s(multi_beam_idx), ...
            r01_p(multi_beam_idx), r12_p(multi_beam_idx), ...
            exp_term(multi_beam_idx));
    end
    
    % 对不需要多光束模型的波长点使用双光束计算
    double_beam_idx = ~use_multi_beam;
    if any(double_beam_idx)
        R(double_beam_idx) = cal_R_db(...
            r01_s(double_beam_idx), r12_s(double_beam_idx), ...
            r01_p(double_beam_idx), r12_p(double_beam_idx), ...
            exp_term(double_beam_idx));
    end
    
    % 输出使用情况统计（调试用）
    if false  % 设为true以显示调试信息
        fprintf('自适应模型使用情况: 多光束%.1f%%, 双光束%.1f%%\n', ...
            100*sum(multi_beam_idx)/length(waveNum), ...
            100*sum(double_beam_idx)/length(waveNum));
    end
end