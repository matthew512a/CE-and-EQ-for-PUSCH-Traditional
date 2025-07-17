clc; clear; close all;

%% Thông số
n_Subcarriers = 3276;       % số lượng sóng mang con (subcarriers)
data_Symbols = 13;          % Số lượng 
pilot_Symbols = 1;          % Số lượng
bps = 4;                    % bit per symbols (data Symbols) (16QAM)
M = 2^bps;                  
u = 1;                      % numerology
SCS = 30000;                % Subcarier Spacing
K = 4;
Nfft = 2048 * K * 2.^(-u);  % Số mẫu IFFT theo chuẩn 3GPP (4096)
Ncp0 = 352;                 % tại l = 0
Ncp = 288;                  % tại l = 1:13

%% Tạo bit ngẫu nhiên và reshape
bits = randi([0 1], n_Subcarriers, data_Symbols, bps);     % size: [3276 x 13 x 4]
disp(bits);
bits_reshaped = reshape(bits, [], bps);                 % size: [3276*13 x 4]
disp(bits_reshaped);

%% Điều chế 16QAM - 3GPP
modulated_symbols = qam16Modulate3gpp(bits_reshaped);  % vector [3276*13 x 1]
disp(modulated_symbols);

%% Bước 3: Tạo lưới tài nguyên và chèn DMRS
pusch_grid = zeros(n_Subcarriers, 14);   % complex grid
%disp(pusch_grid)
DMRS_Signal = readmatrix("dmrs_signal.xlsx", 'Range', 'A1:A3276');
%size(DMRS_Signal)
dmrs_col = 4;  % MATLAB indexing starts at 1
pusch_grid(:, dmrs_col) = DMRS_Signal;
data_cols = setdiff(1:14, dmrs_col);    % Loại bỏ cột 4
k = 1;

for col = data_cols
    pusch_grid(:, col) = modulated_symbols(k:k + n_Subcarriers - 1);
    k = k + n_Subcarriers;
end
%fprintf("Đã gán xong vào grid. Kích thước: %d x %d\n", size(pusch_grid));
%disp(pusch_grid(1:10, :))

%% Bước 4: OFDM modulation + thêm CP
% Thực hiện IFFT cho từng OFDM symbol, chèn Cyclic Prefix (CP)
% TODO: ifft(), chọn CP cuối mỗi symbol (cuối 144 mẫu), ghép lại chuỗi tín hiệu
tx_signal = [];                             % Tín hiệu sau OFDM + CP
n_symbol = data_Symbols + pilot_Symbols;    % 14 OFDM Symbol
% Zero padding để map vào Nfft
null_left  = floor((Nfft - n_Subcarriers)/2);
null_right = Nfft - null_left - n_Subcarriers;
for l=0:n_symbol-1
    % 1. Map lưới tần số vào NFFT
    freq_symbol = pusch_grid(:, l+1); % MATLAB index từ 1
    freq_padded = [zeros(null_left,1); freq_symbol; zeros(null_right,1)];
    
    % 2. IFFT
    time_symbol = ifft(ifftshift(freq_padded), Nfft);
    
    % 3. CP cho từng symbol
    if l == 0
        cp_len = Ncp0;
    else
        cp_len = Ncp;
    end
    cp = time_symbol(end - cp_len + 1:end);

    % 4. Ghép CP + OFDM symbol
    tx_symbol = [cp; time_symbol];
    tx_signal = [tx_signal; tx_symbol];
end

% Kiểm tra
fprintf("Tín hiệu phát cuối cùng có %d mẫu\n", length(tx_signal));


%% Bước 5: Mô phỏng kênh truyền
profile = 'TDL-B';         % Hồ sơ trễ chuẩn TDL-B
fd = 70;                   % Doppler shift ước lượng cho khoảng 25m, tốc độ ~30km/h
tm = 100e-9;               % Delay spread: 100ns
Fs = SCS * Nfft;           % SCS = 30 kHz, Nfft = 4096

% Lấy kênh
channel = getChannelInfo(profile, fd, tm, Fs);

% Truyền tín hiệu
rx_signal = channel(tx_signal);
% Thêm nhiễu Gaussian (AWGN) với SNR = 20 dB
% TODO: awgn()
SNR_dB = 20;
rx_signal_awgn = awgn(rx_signal, SNR_dB, 'measured');
noiseVar = 10^(-SNR_dB / 10);           % noise variance (phương sai nhiễu)

%% Bước 6: OFDM demodulation
% Loại bỏ CP, thực hiện FFT để khôi phục tín hiệu miền tần số
% TODO: fft(), tách các symbol từ chuỗi miền thời gian

% Khởi tạo lại lưới tần số sau khi giải điều chế
rx_grid = zeros(n_Subcarriers, n_symbol);  % kích thước giống grid
% Tính độ dài mỗi symbol (bao gồm cả CP)
symbol_lengths = zeros(1, n_symbol);
for l = 0:n_symbol-1
    if l == 0
        symbol_lengths(l+1) = Nfft + Ncp0;
    else
        symbol_lengths(l+1) = Nfft + Ncp;
    end
end
% Chỉ số bắt đầu mỗi symbol
start_idx = cumsum([1, symbol_lengths(1:end-1)]);
% Demodulation
for l = 0:n_symbol-1
    % Lấy symbol có CP
    idx_start = start_idx(l+1);
    idx_end = idx_start + symbol_lengths(l+1) - 1;
    rx_symbol_with_cp = rx_signal_awgn(idx_start:idx_end);

    % Loại bỏ CP
    if l == 0
        cp_len = Ncp0;
    else
        cp_len = Ncp;
    end
    rx_symbol = rx_symbol_with_cp(cp_len+1:end);

    % FFT
    freq_symbol = fftshift(fft(rx_symbol, Nfft));

    % Cắt phần giữa để thu lại 3276 subcarriers
    rx_grid(:, l+1) = freq_symbol(null_left+1 : null_left+n_Subcarriers);
end

%% Bước 7: Giải điều chế hardbit
% Lấy các cột dữ liệu từ grid đã giải điều chế, giải QAM để ra bit
% TODO: qamdemod(), de2bi(), reshape()
rx_data_grid = rx_grid(:, data_cols);       % [3276 × 13]

% Chuẩn bị vector đầu vào cho demodulation
rx_symbols = reshape(rx_data_grid, [], 1);  % [3276×13 x 1] vector

% Hàm giải điều chế hardbit 16QAM theo chuẩn 3GPP Gray mapping
rx_bits = qam16Demodulate3gpp(rx_symbols);  % [3276×13 x 4] ma trận bit

%% Bước 8: Tính MSE và BER
tx_data_symbols = reshape(modulated_symbols, n_Subcarriers, data_Symbols);

% Symbol thu sau OFDM + AWGN
rx_data_symbols = rx_data_grid;  % đã reshape đúng ở bước 7

% Tính MSE
mse = mean(abs(tx_data_symbols(:) - rx_data_symbols(:)).^2);
fprintf("MSE (mean squared error) = %.5f\n", mse);

% Bit đầu vào gốc
tx_bits = bits;

% Bit đầu ra thu được sau giải điều chế
bit_errors = sum(rx_bits(:) ~= tx_bits(:));
total_bits = numel(tx_bits);
ber = bit_errors / total_bits;

fprintf("BER (bit error rate) = %.6f\n", ber);
