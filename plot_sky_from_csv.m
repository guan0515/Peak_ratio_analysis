function plot_sky_from_csv(csv_file, t_local, prn_ratios, show_all)
% plot_sky_from_csv 根據提供的 CSV 檔案畫出衛星的 Sky Plot (天空圖)
%
%   用法:
%       plot_sky_from_csv(csv_file)
%       plot_sky_from_csv(csv_file, t_local)
%       plot_sky_from_csv(csv_file, t_local, prn_ratios)
%       plot_sky_from_csv(csv_file, t_local, prn_ratios, show_all)
%
%   輸入參數:
%       csv_file   - 字串 (String) 或字元陣列 (char array)，指定 CSV 檔案的名稱或路徑。
%                    CSV 必須包含欄位：PRN, Elevation_deg, Azimuth_deg。
%       t_local    - (選填) datetime 物件，用來標示在圖片標題上的當地時間。
%                    若未提供或給定空矩陣 []，標題將不會顯示時間資訊。
%       prn_ratios - (選填) 1xN 或 Nx1 的數值陣列，代表各個 PRN 的平均 Peak Ratio。
%                    陣列的索引 (Index) 對應 PRN 號碼 (例如 prn_ratios(13) 代表 PRN 13)。
%                    如果有提供這個參數，點的顏色就會依照 Ratio 大小來塗色 (Color-coded)。
%       show_all   - (選填) 布林值 (true/false)。
%                    true : 顯示所有衛星 (包含小於 0 度，預設值)。
%                    false: 只顯示仰角大於 0 度的衛星。

    % 檢查輸入參數數量，處理時間標註與 Ratio 的預設行為
    has_time = (nargin >= 2) && ~isempty(t_local);
    has_ratios = (nargin >= 3) && ~isempty(prn_ratios);
    if nargin < 4 || isempty(show_all)
        show_all = true; % 預設顯示全部
    end

    % 檢查指定的 CSV 檔案是否存在，避免讀取錯誤
    if ~isfile(csv_file)
        error('找不到指定的檔案：%s', csv_file);
    end

    % 讀取 CSV 檔案內容轉換成 MATLAB 的 table 格式
    data = readtable(csv_file);

    % ---------------------------------------------------------
    % 1. 篩選資料
    % ---------------------------------------------------------
    if show_all
        data_filtered = data;          % 不過濾，直接使用全部資料
    else
        idx = data.Elevation_deg > 0;
        data_filtered = data(idx, :);  % 只取出仰角大於 0 度的衛星
    end

    % 提取繪圖所需的三個主要變數
    elev_list = data_filtered.Elevation_deg; % 衛星仰角 (Elevation) 列表
    azim_list = data_filtered.Azimuth_deg;   % 衛星方位角 (Azimuth) 列表
    prn_list_str = data_filtered.PRN;        % 衛星編號 (PRN) 列表，此時為字串格式 (例如 'G13')

    % ---------------------------------------------------------
    % 2. 資料前處理：將 PRN 從字串 (例如 'G13') 轉成數值 (例如 13)
    % ---------------------------------------------------------
    num_sats = height(data_filtered);        % 取得過濾後的衛星總數
    prn_num = zeros(num_sats, 1);            % 初始化儲存數值 PRN 的陣列
    
    for i = 1:num_sats
        % 判斷 prn_list_str 的資料型態以正確提取字串
        if iscell(prn_list_str)
            str_val = prn_list_str{i};       % 如果是 cell array
        else
            str_val = char(prn_list_str(i)); % 如果是 string array
        end
        % 去除字首的 'G' (GPS) 並將剩餘字串轉成數值 (double)
        prn_num(i) = str2double(strrep(str_val, 'G', ''));
    end

    % ---------------------------------------------------------
    % 3. 建立並繪製 Sky Plot (天空圖)
    % ---------------------------------------------------------
    % 建立一個新的 Figure 視窗並設定大小與位置
    if show_all
        figure('Position', [100, 100, 800, 700], 'Name', 'Sky Plot (All Elevations)');
    else
        figure('Position', [100, 100, 800, 700], 'Name', 'Sky Plot (> 0 deg)');
    end
    ax = polaraxes; % 建立極座標軸 (Polar Axes) 用於畫 Sky plot
    hold on;        % 保持畫布以允許在同一個圖上疊加多個元素

    % 極座標轉換：
    % 方位角 (Azimuth)：需從度 (degree) 轉換為弧度 (radian)，以符合 polarscatter 的輸入要求
    theta = deg2rad(azim_list);
    % 仰角 (Elevation)：在 Sky plot 中，天頂 (90度仰角) 在中心，地平線 (0度仰角) 在最外圈。
    % 所以半徑 r = 90 - 仰角
    r = 90 - elev_list;

    % ---------------------------------------------------------
    % 3.5 準備顏色陣列 (Color-coded)
    % ---------------------------------------------------------
    if has_ratios
        colors = zeros(num_sats, 1);
        for i = 1:num_sats
            p = prn_num(i);
            % 確保 PRN 號碼合法且不超過 prn_ratios 的長度
            if p > 0 && p <= length(prn_ratios)
                colors(i) = prn_ratios(p);
            else
                colors(i) = NaN; % 如果沒有對應的資料則設為 NaN (留白)
            end
        end
        
        % 畫出所有衛星的圓點，使用 colors 陣列上色
        polarscatter(theta, r, 100, colors, 'filled', 'MarkerEdgeColor', 'k');
        
        % 加入 Colorbar 並設定 Colormap
        cb = colorbar;
        colormap('turbo'); % 使用與 Heatmap 一致的 colormap
        ylabel(cb, 'Average Peak Ratio', 'FontSize', 12, 'FontWeight', 'bold');
        
        % 將顏色的下限固定在 1 (因為 Ratio 最小通常為 1)
        clim([1, max(max(colors), 1.01)]);
    else
        % 沒有提供 Ratio 時，統一畫成藍色
        polarscatter(theta, r, 100, 'b', 'filled', 'MarkerEdgeColor', 'k');
    end

    % ---------------------------------------------------------
    % 4. 加上衛星編號 (PRN) 標籤
    % ---------------------------------------------------------
    for i = 1:num_sats
        % 在每個點的旁邊 (半徑往外推 4 個單位) 加上文字標籤 (例如 'G13')
        % 往外推是為了避免文字和資料點重疊
        text(theta(i), r(i)+10, sprintf('G%02d', prn_num(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    end

    % ---------------------------------------------------------
    % 5. 設定極座標軸 (Sky Plot) 的外觀與刻度
    % ---------------------------------------------------------
    ax.ThetaDir = 'clockwise';        % 角度方向設為順時針 (符合方位角習慣)
    ax.ThetaZeroLocation = 'top';     % 將 0 度 (正北) 放置在畫布的最上方
    
    if show_all
        ax.RLim = [0 180];                % 半徑限制改為 0 到 180 (對應仰角 90 度到 -90 度)
        ax.RTick = [0 30 60 90 120 150 180];
        ax.RTickLabel = {'90^\circ', '60^\circ', '30^\circ', '0^\circ (Horizon)', '-30^\circ', '-60^\circ', '-90^\circ'}; % 顯示正負仰角
        
        % 畫出一條紅線標示地平線 (0度仰角 -> 半徑 90)
        th_circle = linspace(0, 2*pi, 200);
        polarplot(th_circle, repmat(90, 1, 200), 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    else
        ax.RLim = [0 90];                 % 半徑限制為 0 到 90 (中心點是 0 對應仰角 90，最外圈是 90 對應仰角 0)
        ax.RTick = [0 30 60 90];          % 設定半徑方向的刻度環
        ax.RTickLabel = {'90^\circ', '60^\circ', '30^\circ', '0^\circ (Horizon)'}; 
    end

    % ---------------------------------------------------------
    % 6. 設定圖表標題
    % ---------------------------------------------------------
    if show_all
        base_title = 'GPS Sky Plot (All Elevations)';
    else
        base_title = 'GPS Sky Plot (Elevation > 0^\circ)';
    end

    if has_time
        title_str = sprintf('%s\nLocal Time: %s', base_title, datestr(t_local));
    else
        title_str = base_title;
    end
    title(title_str, 'FontSize', 14, 'FontWeight', 'bold');

    hold off; % 釋放畫布
end
