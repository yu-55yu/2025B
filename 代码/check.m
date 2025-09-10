function [use_multi_beam, q_values] = check(thk, n1_complex, n2_complex, waveNum, theta0, epsilon_threshold)
% CHECK_MULTI_BEAM_INTERFERENCE �ж��Ƿ������������
%   �����о� |q| = |r10*r12| * exp(-sigma*2*d) >= epsilon ���ж�
%
% �������:
%   thk: Ĥ�� (��m)
%   n1_complex: ���Ӳ㸴������ (n+ik)
%   n2_complex: �ĵ׸������� (n+ik)
%   waveNum: ���� (cm^-1)
%   theta0: ����� (����)
%   epsilon_threshold: ��ֵ��Ĭ��0.03��ͨ��ȡ[0.01, 0.05]
%
% ���:
%   use_multi_beam: �߼����飬true��ʾ�ò�����Ҫ���Ƕ��������
%   q_values: ��������Ӧ��qֵ�����ڵ��Ժͷ���

    % ����Ĭ����ֵ
    if nargin < 6
        epsilon_threshold = 0.03;  % Ĭ��ȡ�м�ֵ
    end
    
    % ȷ������Ϊ������
    n0 = 1.0;
    n1 = n1_complex(:);
    n2 = n2_complex(:);
    waveNum = waveNum(:);
    
    % ���㲨�� (��m)
    waveLen = 10000 ./ waveNum;
    
    % ��ȡ�������ʵ�ʵ�����鲿
    n1_real = real(n1);
    k1_imag = imag(n1);
    n2_real = real(n2);
    k2_imag = imag(n2);
    
    % Ӧ��Snell���ɼ��������
    sin_theta1 = n0 * sin(theta0) ./ n1;
    cos_theta1 = sqrt(1 - sin_theta1.^2);
    sin_theta2 = n0 * sin(theta0) ./ n2;
    cos_theta2 = sqrt(1 - sin_theta2.^2);
    
    % ����Fresnel����ϵ�� (sƫ���pƫ��)
    r01_s = (n0*cos(theta0) - n1.*cos_theta1) ./ (n0*cos(theta0) + n1.*cos_theta1);
    r12_s = (n1.*cos_theta1 - n2.*cos_theta2) ./ (n1.*cos_theta1 + n2.*cos_theta2);
    r01_p = (n1*cos(theta0) - n0*cos_theta1) ./ (n1*cos(theta0) + n0.*cos_theta1);
    r12_p = (n2.*cos_theta1 - n1.*cos_theta2) ./ (n2.*cos_theta1 + n1.*cos_theta2);
    
    % ����ƽ������ϵ����ģ
    r01_avg = sqrt((abs(r01_s).^2 + abs(r01_p).^2) / 2);
    r12_avg = sqrt((abs(r12_s).^2 + abs(r12_p).^2) / 2);
    
    % ��������ϵ�� �� = 4��k/�� (��λ: ��m^-1)
    alpha = 4 * pi * k1_imag ./ waveLen;
    
    % ����� = 2�� (����б����ʱ��ʵ�ʹ��)
    % �ڱ�Ĥ�е�ʵ�ʴ���������Ҫ���������
    actual_path_factor = real(cos_theta1);  % ʵ�ʹ����������
    sigma = 2 * alpha ./ actual_path_factor;
    
    % ����˥������ exp(-��*2d)
    % 2d��ʾ���ڱ�Ĥ������һ�εľ���
    attenuation_factor = exp(-sigma * 2 * thk);
    
    % �����о�ֵ |q| = |r10*r12| * exp(-��*2d)
    q_values = r01_avg .* r12_avg .* attenuation_factor;



figure; 
plot(waveNum, q_values, 'LineWidth', 1.2);
hold on;


e_r = [0.6350, 0.0780, 0.1840]; % ���� MATLAB ��  ɫ

yline(epsilon_threshold, 'LineStyle', '--', 'Color',e_r, 'LineWidth', 1.2);

hold off; % ������ɣ��ر� hold

xlabel('���� (cm^{-1})');
ylabel('q ֵ');
grid on; 
xlim([min(waveNum), max(waveNum)]);

legend('q ֵ', sprintf('��ֵ ��=%.2f', epsilon_threshold));
    
    % �ж��Ƿ���Ҫ���Ƕ��������
    use_multi_beam = q_values >= epsilon_threshold;
    
    % ���ͳ����Ϣ����ѡ��
    if nargout ~= 0
        fprintf('\n=== ����������жϽ�� ===\n');
        fprintf('��ֵ �� = %.3f\n', epsilon_threshold);
        fprintf('qֵ��Χ: [%.4f, %.4f]\n', min(q_values), max(q_values));
        fprintf('��Ҫ�����ģ�͵Ĳ�����: %d/%d (%.1f%%)\n', ...
            sum(use_multi_beam), length(use_multi_beam), ...
            100*sum(use_multi_beam)/length(use_multi_beam));
        
        % �ҳ��ٽ粨������
        transition_idx = find(diff(use_multi_beam) ~= 0);
        if ~isempty(transition_idx)
            fprintf('ģ��ת�������� (cm^-1): ');
            fprintf('%.1f ', waveNum(transition_idx));
            fprintf('\n');
        end
    end
end

function R = cal_R_mul(r01_s, r12_s, r01_p, r12_p, exp_term)
% CALCULATE_R_MULTI_BEAM ���������ģ�ͣ�������Airy��ʽ��
    % sƫ��
    numerator_s = abs(r01_s + r12_s .* exp_term).^2;
    denominator_s = abs(1 + r01_s .* r12_s .* exp_term).^2;
    R_s = numerator_s ./ denominator_s;
    
    % pƫ��
    numerator_p = abs(r01_p + r12_p .* exp_term).^2;
    denominator_p = abs(1 + r01_p .* r12_p .* exp_term).^2;
    R_p = numerator_p ./ denominator_p;
    
    % ��ƫ���ķ����ʣ�s��p��ƽ����
    R = real((R_s + R_p) / 2);
end

function R = compute_R_adaptive(thk, n1_complex, n2_complex, waveNum, theta0, epsilon_threshold)
% COMPUTE_R_ADAPTIVE ����Ӧѡ�����ģ�͵ķ����ʺ���
%   �������������Զ�ѡ��˫��������������ģ��
%
% ���������
%   thk: Ĥ�� (��m)
%   n1_complex: ���Ӳ㸴������
%   n2_complex: �ĵ׸�������
%   waveNum: ���� (cm^-1)
%   theta0: ����� (����)
%   epsilon_threshold: ������о���ֵ (Ĭ��0.03)

    if nargin < 6
        epsilon_threshold = 0.03;
    end
    
    % �ж��Ƿ���Ҫ�����ģ��
    [use_multi_beam, q_values] = check_multi_beam_interference(thk, n1_complex, n2_complex, waveNum, theta0, epsilon_threshold);
    
    % ׼����������Ĳ���
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
    
    % ��ʼ���������
    R = zeros(size(waveNum));
    
    % ����Ҫ�����ģ�͵Ĳ�����ʹ�ö��������
    multi_beam_idx = use_multi_beam;
    if any(multi_beam_idx)
        R(multi_beam_idx) = cal_R_mul(...
            r01_s(multi_beam_idx), r12_s(multi_beam_idx), ...
            r01_p(multi_beam_idx), r12_p(multi_beam_idx), ...
            exp_term(multi_beam_idx));
    end
    
    % �Բ���Ҫ�����ģ�͵Ĳ�����ʹ��˫��������
    double_beam_idx = ~use_multi_beam;
    if any(double_beam_idx)
        R(double_beam_idx) = cal_R_db(...
            r01_s(double_beam_idx), r12_s(double_beam_idx), ...
            r01_p(double_beam_idx), r12_p(double_beam_idx), ...
            exp_term(double_beam_idx));
    end
    
    % ���ʹ�����ͳ�ƣ������ã�
    if false  % ��Ϊtrue����ʾ������Ϣ
        fprintf('����Ӧģ��ʹ�����: �����%.1f%%, ˫����%.1f%%\n', ...
            100*sum(multi_beam_idx)/length(waveNum), ...
            100*sum(double_beam_idx)/length(waveNum));
    end
end