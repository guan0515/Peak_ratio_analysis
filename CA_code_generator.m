function CA_code = CA_code_generator(format)
%寫一個生成37種C/A code的function
%輸出: 一個37*1023的矩陣
%輸入: 一個G2的多項式係數表

if nargin < 1
    format = '01';   % 預設輸出 0/1
end

%定義長度
m = 10; % 總共10-bit
n = 2^m - 1; % code長度:1023
G2_poly_table = [2 6;3 7;4 8;5 9;1 9;2 10;1 8;2 9;3 10;2 3;3 4;5 6;6 7;7 8;8 9;9 10;1 4;2 5;3 6;4 7;5 8;6 9;1 3;4 6;5 7;6 8;7 9;8 10;1 6;2 7;3 8;4 9;5 10;4 10;1 7;2 8;4 10];
r = size(G2_poly_table,1); % 所有CA code種類數量

%定義暫存器初始值
G1_reg = ones(1,m); % G1暫存器初始值
G2_reg = ones(1,m); % G2暫存器初始值

%定義序列空間
G1_seq = zeros(1,n); % G1 code序列
G2_seq = zeros(1,n); % G2 code序列
G2I_seq = zeros(r,n); % G2I code序列
CA_code_seq = zeros(r,n); % CA code序列


% Gold code生成
for k = 1:r
    for i = 1:n
        G1_seq(i) = G1_reg(end); % G1 輸出當前數值
        G2_seq(i) = G2_reg(end); % G2 輸出當前數值
        
        G1_feedback = mod(G1_reg(3) + G1_reg(10), 2); 
        G2_feedback = mod(G2_reg(2) + G2_reg(3) + G2_reg(6) + G2_reg(8) + G2_reg(9) + G2_reg(10), 2);
        %(0+0)/2 餘0, (1+1)/2 餘0, (1+0)/2 餘1，結果和XOR一樣(相同為0，相異為1)

        G2I_seq(i) = mod(sum(G2_reg(G2_poly_table(k,:))), 2);
        CA_code_seq(k,i) = xor(G1_seq(i),G2I_seq(i));

        G1_reg = [G1_feedback, G1_reg(1:end-1)]; % 加入feedback，形成新數列
        G2_reg = [G2_feedback, G2_reg(1:end-1)]; % 加入feedback，形成新數列
        
        
    end
end

% ===============================
%  Output format conversion
% ===============================
switch upper(format)
    case 'NRZ'
        CA_code = 2 * CA_code_seq - 1;     % 0→-1, 1→+1
    otherwise
        CA_code = CA_code_seq;             % 0/1
end





% % 檢查前 10 位八進位值
% for k = 1:r
% bin_code = num2str(CA_code_seq(k,1:10));
% decimal_value = bin2dec(bin_code);
% octal_value = dec2base(decimal_value, 8);
% 
% disp(['C/A code No.',num2str(k)]);
% disp(['前 10 位 C/A 碼 (二進位): ', bin_code]);
% disp(['轉換後八進位值: ', num2str(octal_value)]);
% fprintf('\n');  % 輸出一個空行
% end
% 
% %生成檔案
% filename = 'CA_code.txt';
% 
% %存入檔案
% fid = fopen(filename, 'w');
% 
% % 寫入標題行
% fprintf(fid, 'CA Code No.\tFirst 10 Bits\tOctal Value\n');
% 
% for k = 1:r
%     bin_code = num2str(CA_code_seq(k,1:10)); % 取前10位二進位
%     decimal_value = bin2dec(bin_code); % 轉換成十進位
%     octal_value = dec2base(decimal_value, 8); % 轉換成八進位
% 
%     % 將CA code種類號碼、前10個bit、八進位數值寫入檔案
%     fprintf(fid, '%d\t%s\t%s\n', k, bin_code, octal_value); %*1
% end
% fclose(fid);
% 
% disp('succeeded')


end

%*1
% '%d': 這個格式符號表示將一個整數轉換為數字格式。在這裡，它會將 k（即 CA code 種類號碼）轉換為整數並寫入檔案。
% '\t': 這是制表符（tab），它會在不同欄位之間加入空白，讓輸出的檔案有更清晰的結構。每個 \t 會插入一個水平制表符，讓每一欄位之間有間隔。
% '%s': 這個格式符號表示將一個字串（string）輸出。這裡，bin_code 和 octal_value 都是字串，因此會使用 %s 格式符號來將它們寫入檔案。
% '\n': 這個格式符號表示換行符號，讓每一行的內容之後換行。