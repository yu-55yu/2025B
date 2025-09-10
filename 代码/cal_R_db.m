function R = cal_R_db(r01_s, r12_s, r01_p, r12_p, exp_term)
% CALCULATE_R_DOUBLE_BEAM 双光束干涉模型
    r_s_total = r01_s + r12_s .* exp_term;
    r_p_total = r01_p + r12_p .* exp_term;
    R_s = abs(r_s_total).^2;
    R_p = abs(r_p_total).^2;
    R = real((R_s + R_p) / 2);
end