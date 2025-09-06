clc;
clear;
close all;

% --- 文件和测量配置
file1.path = '附件1.xlsx';
file1.angle = 10; % 入射角 (度)
file2.path = '附件2.xlsx';
file2.angle = 15; % 入射角 (度)

% --- 拟合全局配置
config.waveN_fit_min = 1200; % 拟合起始波数 (cm??)
config.n1_init = 2.58;       % 外延层折射率初值 (用于FFT估算厚度)
config.n2_real_init = 2.55;  % 衬底折射率初值
config.k1_A_init = 0.001;    % 外延层消光系数参数A初值
config.k1_B_init = 2.0;      % 外延层消光系数参数B初值
% Sellmeier模型 B1, B2, B3, C1, C2, C3 参数初值
config.selParam_init = [5.5, 0.2, 0.05, 0.027, 100, 0.01]; 

% --- 步骤 1: 使用FFT估算厚度初值
disp('=====================================================');
disp('步骤 1: 计算厚度初值');
data1 = readmatrix(file1.path);
thk_init1 = fft_thk_estimate(data1(:,1), data1(:,2)/100, config.n1_init, config.waveN_fit_min, file1.angle);

data2 = readmatrix(file2.path);
thk_init2 = fft_thk_estimate(data2(:,1), data2(:,2)/100, config.n1_init, config.waveN_fit_min, file2.angle);

% 使用两个角度估算结果的平均值作为全局初值
config.thk_init = mean([thk_init1, thk_init2]);
fprintf('文件1 (10°) FFT估算厚度: %.2f μm\n', thk_init1);
fprintf('文件2 (15°) FFT估算厚度: %.2f μm\n', thk_init2);
fprintf('<strong>平均厚度初值: %.2f μm</strong>\n', config.thk_init);
disp('=====================================================');

% --- 步骤 2: 对每个文件进行独立的全局拟合
disp(' ');
disp('步骤 2: 对每个文件进行拟合');
disp('--> 开始处理附件1 (10°)...');
[result1] = process(data1, file1.angle, config, file1.path);
disp('--> 开始处理附件2 (15°)...');
[result2] = process(data2, file2.angle, config, file2.path);
disp('=====================================================');


% --- 步骤 3: 分析和比较两个角度的拟合结果
disp(' ');
disp('步骤 3: 综合分析结果');
analyze_res(result1, result2, file1.angle, file2.angle);
disp('=====================================================');




%% 数据处理和拟合主函数
function [result] = process(data, angle, config, output_filepath)
    waveNum = data(:, 1);       % 波数 (cm??)
    R = data(:, 2) / 100;       % 反射率 (0-1)
    waveLen_full = 10000 ./ waveNum; % 波长 (μm)

    % 根据设定的波数范围筛选用于拟合的数据
    filter = waveNum >= config.waveN_fit_min;
    waveNum_fit = waveNum(filter);
    R_fit = R(filter);
    
    % 参数顺序: [厚度, 衬底n2, Sellmeier(6个), k参数(2个)]
    x0 = [config.thk_init, config.n2_real_init, config.selParam_init, config.k1_A_init, config.k1_B_init];
    lb = [config.thk_init*0.8, 2.0, 0.1, 0.001, 0.0001, 0.0001, 0.1, 0.001, 0, 0];
    ub = [config.thk_init*1.2, 3.5, 20, 10, 5, 2, 150, 20, 0.01, 4];

    model_type = 'real_substrate'; 
    [x_optimal, R_squared_fit] = global_fit(x0, lb, ub, waveNum_fit, R_fit, angle, model_type);

    result.thk = x_optimal(1);
    result.n2 = x_optimal(2);
    result.selParam = x_optimal(3:8);
    result.k1Param = x_optimal(9:10);
    result.n1_complex_full = cal_n_sellmeier(waveLen_full, result.selParam, result.k1Param); % <--- 已修正函数名

    fprintf('    最终厚度: %.2f μm\n', result.thk);
    fprintf('    最终外延层 n1 (在 6μm): %.3f + %.4fi\n', real(interp1(waveLen_full, result.n1_complex_full, 6)), imag(interp1(waveLen_full, result.n1_complex_full, 6)));
    fprintf('    最终衬底 n2: %.3f\n', result.n2);
    fprintf('    在拟合区域的拟合优度 R?: %.4f\n', R_squared_fit);

    theta0_rad = angle * pi / 180;
    R_fit_full = compute_R(result.thk, result.n1_complex_full, result.n2, waveNum, theta0_rad, 1);
    
    try
        R_fit_per = R_fit_full * 100;
        header = {'拟合反射率'};
        writecell(header, output_filepath, 'Sheet', 1, 'Range', 'C1');
        writematrix(R_fit_per, output_filepath, 'Sheet', 1, 'Range', 'C2');
        fprintf('    成功将拟合结果保存到: %s\n', output_filepath);
    catch ME
        fprintf('    保存Excel文件时出错: %s\n', ME.message);
    end

    plot_fit_res(waveNum, R, R_fit_full, angle, config.waveN_fit_min);
    plot_optical_constants(waveLen_full, waveNum, result.n1_complex_full, '外延层 (Sellmeier)');
end


