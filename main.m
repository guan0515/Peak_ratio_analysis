%% ====================================
%       Main section
%  ====================================
%   函數執行區域
%   目標: 
%   1. 執行acquisition
%   2. 畫圖demo 結果 (20260429)
%   3. 觀察後續ratio變化
%
%   label fontsize = 14
%   title / legend fontsize = 15
%   linewidth = 1.1
%   plot, title, label, legend, grid on

%% ====================================
%       clc; clear; close all;
%  ====================================
    clc; clear; close all;


%% ====================================
%       Load data
%  ====================================
    iq_file = 'signal_source_20260313.mat';
    load(iq_file);
    
%% ====================================
%       Remove DC offset 
%  ====================================
    IQ_noDC = IQ - mean(IQ);
    
%% ====================================
%       Parameter setting
%  ====================================
    fs           = 2e6;    % sampling rate = 2 MHz
    doppler_step = 500;    % step = 500 Hz
    N_block      = 10;     % 想要多長段數的訊號
    coh_block    = 1;      % 決定acquisition策略

% %% ========================================================================
% %       Acquisition block (Non-coherent: 10 ms) 
% %  ========================================================================
%     blockId      = 1;      % 決定從整串IQ中的"甚麼位置"開始向後取N_block長資料
%     acq_results  = fft_acquisition(IQ_noDC, fs, doppler_step, blockId, N_block, coh_block, iq_file);
%     disp(acq_results)

% %% ========================================================================
% %       Acquisition block (Non-coherent: 10 ms) 
% %  ========================================================================
%     for i = 1:10:101
%         blockId      = i;      % 決定從整串IQ中的"甚麼位置"開始向後取N_block長資料
%         acq_results  = fft_acquisition(IQ_noDC, fs, doppler_step, blockId, N_block, coh_block, iq_file);
%         disp(acq_results)
%     end


%% ========================================================================
%       Analysis the ratio vs time (ms) 
%  ========================================================================
    time_ms = [];
    all_ratios = [];
    set = 10; % 檢查組數(總共做幾次acquisition)
    for i = 1:N_block:N_block * (set-1) + 1
        blockId      = i;      % 決定從整串IQ中的"甚麼位置"開始向後取N_block長資料
        acq_results  = fft_acquisition(IQ_noDC, fs, doppler_step, blockId, N_block, coh_block, iq_file);
        time_ms(end+1) = (blockId - 1);
        all_ratios(end+1,:) = acq_results.Ratio';
    end
%% ========================================================================
%       Plot the ratio vs time (ms) 
%  ========================================================================
    %detected_mask = any(all_ratios > 1.35, 1);
    %PRN_to_plot = find(detected_mask);

    %PRN_to_plot = [5, 6, 7, 11, 12, 13, 15, 19, 20, 21, 25, 29, 30]; % 仰角大於等於0
    %PRN_to_plot = [5, 6, 11, 13, 15, 20, 21, 29, 30]; % 仰角大於等於15
    %PRN_to_plot = [7, 12, 19, 25]; % 仰角大於等於0 & 小於15

    % exclude_prns = [5, 6, 7, 11, 12, 13, 15, 19, 20, 21, 25, 29, 30];
    % PRN_to_plot  = setdiff(1:37, exclude_prns);
    % 
    % % --- plot
    % figure('Position',[100,100,1200,600]);
    % hold on;
    % colors = turbo(length(PRN_to_plot));
    % for k = 1:length(PRN_to_plot)
    %     prn = PRN_to_plot(k);
    %     plot(time_ms, all_ratios(:, prn), '-o', 'LineWidth', 1.1, 'DisplayName', sprintf('PRN %d', prn), 'Color', colors(k, :));
    % end
    % 
    % yline(1.35, '--k', 'Threshold = 1.35','DisplayName', sprintf('Threshold'), 'LineWidth', 1.5, 'LabelHorizontalAlignment','left');
    % xlabel('Times (ms)', 'FontSize', 14);
    % ylabel('Ratio (Peak / 2nd Peak)', 'FontSize', 14);
    % title('PRN Ratio Trend vs Time', 'FontSize', 15);
    % legend('Location','eastoutside', 'FontSize', 10);
    % grid on;

%% ========================================================================
%       Heatmap the ratio vs time (ms) 
% %  ========================================================================  
%     figure('Position',[100,100,1200,600]);
%     imagesc(time_ms, 1:37, all_ratios');
%     colormap('turbo');
%     cb = colorbar;
%     ylabel(cb, 'Ratio', 'FontSize', 13)
%     %clim([1 max(all_ratios(:))*1.05]);
%     clim([1, prctile(all_ratios(:), 95)]);
%     hold on;
%     yline((1:37)-0.5, 'w-', 'LineWidth', 0.3)
% 
%     xlabel('Times (ms)', 'FontSize', 14);
%     ylabel('PRN', 'FontSize', 14);
%     title('PRN Ratio Heatmap vs Time', 'FontSize', 15);
%     yticks(1:37);
%     yticklabels(1:37);
    
%% ========================================================================
%       Sky plot
%  ======================================================================== 
    t_local = datetime(2026,3,13,19,46,0,'TimeZone','Asia/Taipei');
    
    % 計算每個 PRN 的平均 Ratio (針對時間維度平均)
    mean_ratios = mean(all_ratios, 1);
    
    % 畫出天空圖，並根據平均 Ratio 上色
    % 第四個參數設定是否畫出所有衛星 (true: 畫出全部含小於0度, false: 只畫大於0度)
    show_all_elevations = true;
    plot_sky_from_csv('20260313satellite_positions.csv', t_local, mean_ratios, show_all_elevations);
    
    % 加上副標題顯示檢查組數與總時間
    % 1個 N_block 代表 1 ms，所以總時間 = set * N_block / 1000 秒
    total_sec = (set * N_block) / 1000;
    sub_str = sprintf('Average of %d times acquisition (%g sec)', set, total_sec);
    subtitle(sub_str, 'FontSize', 12, 'FontWeight', 'bold');