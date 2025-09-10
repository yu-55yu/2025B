function R = cal_R_mul(r01_s, r12_s, r01_p, r12_p, exp_term)
% CALCULATE_R_MULTI_BEAM 多光束干涉模型
    r_s_total = (r01_s + r12_s .* exp_term) ./ (1 + r01_s .* r12_s .* exp_term);
    r_p_total = (r01_p + r12_p .* exp_term) ./ (1 + r01_p .* r12_p .* exp_term);
    R_s = abs(r_s_total).^2;
    R_p = abs(r_p_total).^2;
    R = real((R_s + R_p) / 2);
end