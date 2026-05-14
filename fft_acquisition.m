%% ====================================
%       FFT_method_acquisition (function)
%  ====================================
%   FFT頻域法解算用GNSS-SDR接收下來的IQ訊號
%   1. IQ raw data (N*1 format)
%   2. fs: sampling frequency 
%   3. doppler step: 決定doppler軸掃描精度
%   4. blockId: 從第幾段開始取資料
%   5. N_block: 取多長資料 (ms 為單位)
%   6. coh_block: 同相積分長度


% ver4.0 (update: 20260407)
function acq_results = fft_acquisition(IQ, fs, doppler_step, blockId, N_block, coh_block, filename_in)

%% ========================================================================
%       data preprocessing
%  ========================================================================
    [~, name, ~] = fileparts(filename_in); % name

%% ====================================
%       extract ?ms for acquisition
%  ====================================
    Tp        = 1e-3;                      % 1 ms 週期
    Ns_1ms    = round(fs * Tp);            % 每 Tp 有多少 samples
    Ns_ms     = round(fs * Tp * N_block);  % 所有 Tp 有多少 samples
    start_idx = (blockId-1)*Ns_1ms + 1;    % 起始 sample index
    stop_idx  = start_idx + Ns_ms - 1;     % 結束 sample index
    block_IQ  = IQ(start_idx:stop_idx);  % 擷取的訊號
    block_IQ  = block_IQ(:);
    fprintf('this data is No.%d~%d ms\n', blockId-1, blockId-1 + N_block) % 第1個block對應0 ~ 1 ms

%% ====================================
%       整理成每1 ms 一組 總共10組
%  ====================================
    noncoh_block = N_block / coh_block;
    Ns_coh = coh_block * Ns_1ms;
    %block_IQ_group = zeros(noncoh_block, Ns_coh);
    fprintf('Coherent length: %d ms / Noncoherent times: %d\n', coh_block, noncoh_block)

    block_IQ_group = reshape(block_IQ, Ns_coh, noncoh_block);   

    % for g = 1:noncoh_block
    %     start_idx_group = (g-1)*Ns_coh + 1;  % 起始 sample index
    %     stop_idx_gruop  = start_idx_group + Ns_coh - 1;   % 結束 sample index
    %     block_IQ_group(g,:) = block_IQ(start_idx_group:stop_idx_gruop);
    % end
    % block_IQ_group = block_IQ_group.';


%% ========================================================================
%       data acquisition 非同相積分迴圈 分10組
%  ========================================================================

%% ====================================
%       Create result table
%  ====================================
    acq_results = table( ...
        'Size',[0 7], ...
        'VariableNames',{'PRN','StepHz','DelayEst','DopplerEst','Peak (power)','Detected_PRN', 'Ratio'}, ...
        'VariableTypes',{'double','double','double','double','double','double','double'});

%% ====================================
%       Doppler shift removal (Vectorization)
%  ====================================
    doppler_bins = -5000 : doppler_step : 5000 ;   % range = 5000 ~ -5000 Hz
    num_bins     = length(doppler_bins);

    block_mix_group = zeros(Ns_coh, num_bins, noncoh_block);

    for g_idx = 1 : noncoh_block 

        block_IQ_group_index = block_IQ_group(:,g_idx);
        N                    = length(block_IQ_group_index);  % time axis
        t                    = (0:N-1) / fs;   

        doppler_phasors                = exp(-1j * 2*pi * (t.' * doppler_bins));
        block_mix_group(:, :, g_idx) = block_IQ_group_index .* doppler_phasors; % 2000 * 21 (由左至右 -5000~5000)
    end

%% ====================================
%       C/A code generation
%  ====================================
    CA_code = CA_code_generator('NRZ');

%% ====================================
%       fractional code sampling
%  ====================================
    chip_rate = 1.023e6;     % f_c = 1.023MHz
    Nc        = 1023;        % C/A code chips per ms 
    chip_idx  = floor(t * chip_rate);

    % ===== 討論 =====
    % 1. 第n個sample發生在t_n = n / fs
    % 2. 每個chip時間寬度T_c = 1 / fc
    % 3. t / T_c = t * f_c
    % 4. chip(n) = n * f_c / f_s

    chip_idx   = mod(chip_idx, Nc) + 1;
    CA_code_ms = CA_code(:,chip_idx);  % C/A code (fractional code sampling)

    PRN_count  = size(CA_code_ms, 1);
    L          = size(CA_code_ms, 2);

%% ====================================
%       開始計時 (tic)
%  ====================================
    tic;

%% ====================================
%       FFT of received signal
%  ====================================
    FFT_block_mix_group = zeros(L,num_bins,noncoh_block);
    for g_idx = 1:noncoh_block
        FFT_block_mix_group_selected = block_mix_group(:,:,g_idx);
        FFT_block_mix_group(:,:,g_idx) = fft(FFT_block_mix_group_selected, [], 1); % column-wise 
    end
%% ====================================
%       FFT of replica (1 ~ 37)
%  ==================================== 
    FFT_CA_code_ms = fft(CA_code_ms, [], 2); % row-wise

%% ====================================
%       Cross-correlation  (非同相積分)
%  ====================================
    %corr_map_final = zeros(L, num_bins, PRN_count); % 建立 2000×21×37（所有PRN的search plane）
    corr_map_group = zeros(L, num_bins, PRN_count, noncoh_block); % 建立 2000×21×37×N_blocks（所有PRN的search plane）
    for g_idx = 1 : noncoh_block
        block_mix_group_index = FFT_block_mix_group(:, :, g_idx); 

        for prn = 1:PRN_count

            fft_replica_selected = double(FFT_CA_code_ms(prn,:));   % 1×2000 (seq)

            for k = 1:num_bins

                r = block_mix_group_index(:,k); 

                % 每個 delay
                corr_freq = conj(fft_replica_selected.') .* r;
                corr_time = ifft(corr_freq);
                corr_map_group(:,k,prn,g_idx) = abs(corr_time / L).^2; % correlation value (.^2平方 / 單位:power) % 正規化
            end
        end
    end
    
    corr_map_final = mean(corr_map_group,4);
   
%% ====================================
%       結束計時 (toc)  
%  ====================================
    total_fft_time = toc;
    fprintf('FFT-based Acquisition Time: %.4f seconds\n', total_fft_time);
%% ====================================
%       save corr_map  
%  ====================================
    %save([name, '_corr_map.mat'], 'corr_map_group', 'corr_map_final', '-v7.3');


%% ====================================
%       Find peaks  
%  ====================================
    peak_vals       = zeros(PRN_count,1);
    peak_delay_i    = zeros(PRN_count,1);
    peak_doppler_i  = zeros(PRN_count,1);
    ratio_i         = zeros(PRN_count,1);
    SN_ratio_i      = zeros(PRN_count,1);
    CN0_ratio_i     = zeros(PRN_count,1);
    boundary_1ms    = 1:Ns_1ms;
    for prn = 1:PRN_count

        % 單一 PRN 的 correlation plane (2000 × 21)
        corr_plane = corr_map_final(boundary_1ms,:,prn);

        % 最大值 & index
        [peak_vals(prn), idx]    = max(corr_plane(:));
        [delay_idx, doppler_idx] = ind2sub(size(corr_plane), idx);

        peak_delay_i(prn)   = delay_idx - 1;         % delay (samples)
        peak_doppler_i(prn) = doppler_bins(doppler_idx); % doppler (Hz)


        % ---- noise floor estimation (exclude a guard window around peak) ----
        guard_delay   = round(fs / chip_rate);   % 你可以調大/調小
        guard_doppler = (1/Tp)/doppler_step;    % doppler bin guard / 1ms 對應到sinc是1000Hz 一個null

        mask = true(size(corr_plane));
        exclude_idx = mod((delay_idx - guard_delay : delay_idx + guard_delay) - 1, size(corr_plane,1)) + 1;
        f1   = max(1, doppler_idx-guard_doppler); 
        f2   = min(size(corr_plane,2), doppler_idx+guard_doppler);
        mask(exclude_idx, f1:f2) = false;

        noise_samples = corr_plane(mask);                % 排除peak周遭區域值剩下的noise值

        % 找出第二個peak 算第一高點&第二高點的比值
        second_peak = max(noise_samples(:));
        ratio_i(prn) = peak_vals(prn)/second_peak;

        %計算C/N0 (粗估)
        % SN_ratio_i(prn)  = 10*log10((peak_vals(prn) / mean(noise_samples)));
        % Tint = coh_block * 1e-3;
        % CN0_ratio_i(prn) = SN_ratio_i(prn) - 10*log10(Tint);
        % %不適用因為有做non-coherent (公式推算是給1 ms coherent用的)

        %計算RD threshold
        % 不適用，因為涉及Non coherent 沒有close form 去計算Pfa的反推threshold

        detected_PRN = NaN;
        if (ratio_i(prn) > 1.35) % Monte-Carlo法模擬結果
            detected_PRN = prn;
        end

        Trow = table(prn, doppler_step, peak_delay_i(prn), peak_doppler_i(prn), peak_vals(prn), detected_PRN, ratio_i(prn),...
            'VariableNames', {'PRN','StepHz','DelayEst','DopplerEst','Peak (power)','Detected_PRN','Ratio'});

        acq_results = [acq_results; Trow];
    end

% %% ====================================
% %       3D Correlation Plot   
% %  ====================================
%     valid_prns = acq_results.Detected_PRN(~isnan(acq_results.Detected_PRN));
%     %valid_prns = 1:size(corr_map_final,3); %所有PRN
%     for i = 1:length(valid_prns) 
%         prn = valid_prns(i);                            % 指定 PRN 編號
%         corr_plane = corr_map_final(boundary_1ms,:,prn); 
%         delay_axis   = 0:Ns_1ms-1;
% 
%         [Delay, Doppler] = ndgrid(delay_axis, doppler_bins);
%         peak_delay = peak_delay_i(prn);
%         peak_doppler = peak_doppler_i(prn);
% 
%         figure('Position',[100,100,1000,750]);
%         mesh(Delay, Doppler, corr_plane);
%         xlabel('Delay (sample)');
%         ylabel('Doppler (Hz)');
%         zlabel('Correlation Power');
%         title(sprintf('PRN %d | Delay = %d | Doppler = %d Hz | Step = %d Hz', prn, peak_delay, peak_doppler, doppler_step));
%         grid on;
% 
%         cb = colorbar('southoutside');
%         ylabel(cb, 'Power Intensity', 'FontSize',12);
%         colormap("jet");
%         clim([0 max(peak_vals)]);
% 
%         view([-45 35]);
%     end


% %% ====================================
% %       3D Correlation Plot (1~37 GIF) 
% %  ====================================
%     filename = [name, '_correlation.gif'];
% 
%     gif_delay = 0.8;                         % 每一幀的延遲時間 (秒)
% 
%     % figure setting
%     h_fig = figure('Position',[100,100,1000,750],'Visible','on');
%     delay_axis   = 0:Ns_1ms-1;
%     [Delay, Doppler] = ndgrid(delay_axis, doppler_bins);
% 
%     % for loop (1~37)
%     for prn = 1:PRN_count
%         clf(h_fig);
% 
%         corr_plane = corr_map_final(boundary_1ms,:,prn);
%         mesh(Delay, Doppler, corr_plane);
%         hold on; 
%         grid on;
% 
%         cb = colorbar('southoutside'); % 放在下方
%         ylabel(cb, 'Power Intensity', 'FontSize',12);
%         colormap("jet");
%         clim([0 max(peak_vals)]);
% 
%         view([-45 35]); % 設定固定視角
%         xlabel('Delay (sample)','FontSize',14);
%         ylabel('Doppler (Hz)','FontSize',14);
%         zlabel('Correlation Power','FontSize',14);
%         zlim([0 1.1*(max(peak_vals))])
%         peak_delay = peak_delay_i(prn);
%         peak_doppler = peak_doppler_i(prn);
%         title(sprintf('PRN %d | Delay = %d | Doppler = %d Hz | Step = %d Hz', prn, peak_delay, peak_doppler, doppler_step),'FontSize',15);
% 
%         drawnow;
% 
%         % Capture frame
%         frame = getframe(h_fig);
%         img   = frame2im(frame);
%         [imind,cm] = rgb2ind(img,256); % 轉為索引圖片 (GIF 格式)
% 
%         % Write GIF
%         if prn == 1
%             % 第一幀：建立 GIF 檔案
%             imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',gif_delay);
%         else
%             % 後續幀：附加到 GIF 檔案
%             imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',gif_delay);
%         end
%     end
%     close(h_fig);
%     fprintf('GIF 製作完成！檔案已存為: %s\n', filename);

end




