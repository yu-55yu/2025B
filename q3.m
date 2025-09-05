% analyze_epilayer_thickness.m
% 说明：运行前把 附件3.xlsx (10deg) 和 附件4.xlsx (15deg) 放同一目录
% 使用方法：结果会打印出来并绘图。若希望我代跑，请把文件上传。

function analyze_epilayer_thickness()
    files = {'附件3.xlsx','附件4.xlsx'}; % 10deg, 15deg
    angles_deg = [10,15];
    for i=1:2
        fprintf('--- Processing %s (incidence %g°) ---\n', files{i}, angles_deg(i));
        T = readmatrix(files{i});
        wn = T(:,1); % wave number, cm^-1
        Rperc = T(:,2); % percent
        R = Rperc/100;
        % ---- Preprocess: select usable wn range (optional) ----
        % remove NaN and sort
        valid = ~isnan(wn) & ~isnan(R);
        wn = wn(valid); R = R(valid);
        [wn, idx] = sort(wn); R = R(idx);
        % Interpolate to uniform spacing in wn for FFT
        Npoints = 2^nextpow2(length(wn));
        wn_uniform = linspace(min(wn), max(wn), Npoints);
        R_interp = interp1(wn, R, wn_uniform, 'pchip');
        % Detrend and window
        R_detr = detrend(R_interp, 2); % remove slow baseline (2nd order)
        w = hann(length(R_detr))';
        R_win = R_detr .* w;
        % FFT
        Y = fft(R_win);
        P2 = abs(Y/length(Y));
        P1 = P2(1:floor(end/2));
        freq = (0:length(P1)-1)./ (wn_uniform(2)-wn_uniform(1)) / length(wn_uniform); 
        % freq units: cycles per cm^-1 (approx)
        % find peak excluding DC (index 1)
        [pks, locs] = findpeaks(P1, 'SortStr', 'descend');
        if isempty(locs)
            fprintf('  No clear FFT peak found.\n');
            continue;
        end
        main_loc = locs(1);
        f_peak = freq(main_loc);
        % Optical thickness and d (requires n1 and theta1)
        theta0 = angles_deg(i)*pi/180;
        % --- User must supply n1 estimate or model. Here we assume n1 ~ 3.5 as placeholder ---
        n1_guess = 3.45; 
        % compute theta1 via Snell
        n0 = 1;
        theta1 = asin(n0*sin(theta0)/n1_guess);
        OT = f_peak/2; % optical thickness n1*d*cos(theta1)
        d_est = OT / (n1_guess * cos(theta1)); % units: cm^-1? check units: wn in cm^-1 so f_peak is cycles per cm^-1 -> OT in cycles per cm^-1 ...
        % Convert units: since wn in cm^-1, f_peak has unit cycles/(cm^-1) -> OT in (cycles/(cm^-1))/2 -> dimension cm ?
        % Simpler: from formula d (cm) = f_peak/(2*n1*cos(theta1)) * (1)^(?)  -> actually derived earlier units consistent if wn in cm^-1
        % d_est is in cm, convert to um
        d_um = d_est * 1e4; % cm -> um
        fprintf('  FFT peak freq f=%.6g (cycles per cm^-1) -> OT = %.6g (cm^-1 units), estimated d ≈ %.3f um (using n1=%.3g)\n', f_peak, OT, d_um, n1_guess);
        % Simple multi-beam heuristic: peak width and height
        peak_width = estimate_peak_width(P1, main_loc);
        fprintf('  FFT peak width (approx) = %.6g (smaller => higher finesse => more multi-beam)\n', peak_width);
        % Plot
        figure;
        subplot(2,1,1);
        plot(wn, Rperc); xlabel('wave number (cm^{-1})'); ylabel('Reflectance (%)'); title(sprintf('%s raw', files{i}));
        subplot(2,1,2);
        plot(freq, P1); hold on; plot(freq(main_loc), P1(main_loc),'ro'); xlabel('freq (cycles per cm^{-1})'); ylabel('FFT mag'); title('FFT of detrended reflectance');
    end
end

function w = estimate_peak_width(P1, idx)
    % crude: find half-maximum width in index units
    peak = P1(idx);
    half = peak/2;
    left = find(P1(1:idx) < half, 1, 'last');
    if isempty(left), left = 1; end
    right = idx - 1 + find(P1(idx:end) < half, 1, 'first');
    if isempty(right), right = length(P1); end
    w = right - left;
end
