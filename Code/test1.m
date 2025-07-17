    clc; clear; close all;

%% Thông số cấu hình 5G
n_Subcarriers = 3276;       % Number of Subcarriers
nbOFDMforData = 13;         % Number of OFDM Symbol for Data
nbOFDMforDmrs = 1;          % Number of OFDM Symbol for Pilot
posOFDMforDmrs = 4; posOFDMforData = [1 2 3 5 6 7 8 9 10 11 12 13 14];
bps = 4; moduType = '16QAM';
u = 1; SCS = 30e3; K = 4;
Nfft = 2048 * K * 2.^(-u);  % Số mẫu IFFT theo chuẩn 3GPP (4096)
Ncp0 = 352;                 % Number of CP symbols at l = 0
Ncp = 288;                  % Number of CP symbols at l = 1:13


%% Tạo chuỗi bit ngẫu nhiên
bits = randi([0 1], n_Subcarriers*nbOFDMforData*bps, 1);   % [3276*13*4 x 1]
modulated_Syms = nrModuMapper(bits, moduType);
%disp(modulated_Syms(1:10));
% Dùng hàm có sẵn
%bits_matrix = reshape(bits, [], bps);   % [N x 4] matrix
%modulated_Syms = qammod(bits_matrix, 16, 'InputType', 'bit', 'UnitAveragePower', true);
%% Ánh xạ vào lưới 3276 * 14
pusch_grid = zeros(n_Subcarriers, 14);   % complex grid
%disp(pusch_grid);
DMRS_Signal = readmatrix("dmrs_signal.xlsx", 'Range', 'A1:A3276');
%size(DMRS_Signal);
  
pusch_grid(:, posOFDMforDmrs) = DMRS_Signal;

k = 1;
for pos = posOFDMforData
    pusch_grid(:, pos) = modulated_Syms(k:k + n_Subcarriers - 1);
    k = k + n_Subcarriers;
end
fprintf("Đã gán xong vào grid. Kích thước: %d x %d\n", size(pusch_grid));
%disp(pusch_grid(1:10, :));

%% Bước 4: OFDM modulation + thêm CP
% Thực hiện IFFT cho từng OFDM symbol, chèn Cyclic Prefix (CP)
% TODO: ifft(), chọn CP cuối mỗi symbol (cuối 144 mẫu), ghép lại chuỗi tín hiệu
tx_signal = [];                               % Tín hiệu sau OFDM + CP
NSlotSymb = nbOFDMforData + nbOFDMforDmrs;    % Number of Symbol for 1 Slot
% Zero padding để map vào Nfft
null_left  = floor((Nfft - n_Subcarriers)/2);
null_right = Nfft - null_left - n_Subcarriers;
for l=0:NSlotSymb-1
    %fprintf('\n===== Symbol %d =====\n', l);
    % 1. Map lưới tần số vào NFFT
    freq_symbol = pusch_grid(:, l+1); % MATLAB index từ 1
    freq_padded = [zeros(null_left,1); freq_symbol; zeros(null_right,1)];
    
    %fprintf('-> Tín hiệu miền tần số (có padding) - kích thước: %d\n', length(freq_padded));
    %if l < 3
    %    figure; 
    %    plot(real(freq_padded)); title(['Real - Frequency Symbol ', num2str(l)]);
    %    figure; 
    %   plot(imag(freq_padded)); title(['Imag - Frequency Symbol ', num2str(l)]);
    %end

    
    % 2. IFFT
    time_symbol = ifft(ifftshift(freq_padded), Nfft);
    
    %fprintf('-> Sau IFFT (symbol thời gian) - size: %d\n', length(time_symbol));
    %if l < 3
    %    figure; 
    %    plot(real(time_symbol)); title(['Real - Time Symbol ', num2str(l)]);
    %end

    % 3. CP cho từng symbol
    if l == 0
        cp_len = Ncp0;
    else
        cp_len = Ncp;
    end
    cp = time_symbol(end - cp_len + 1:end);
    
    %fprintf('-> CP (length = %d), tổng sau ghép = %d mẫu\n', cp_len, cp_len + length(time_symbol));
    %if l < 3
    %    figure; 
    %    plot(real(cp)); title(['Real - CP Symbol ', num2str(l)]);
    %end

    % 4. Ghép CP + OFDM symbol
    tx_symbol = [cp; time_symbol];
    tx_signal = [tx_signal; tx_symbol];

    % Nếu bạn muốn in vài mẫu đầu:
    %if l < 3
    %    disp('5 giá trị đầu sau ghép CP + Time Symbol:');
    %    disp(tx_symbol(1:5));
    %end
end
    
% Kiểm tra
fprintf("Tín hiệu phát cuối cùng có %d mẫu\n", length(tx_signal));

%% Bước 5: Mô phỏng kênh truyền
profile = 'TDL-B';         % Hồ sơ trễ chuẩn TDL-B
fd = 25;                   % Doppler shift ước lượng cho khoảng 25m, tốc độ ~30km/h
tm = 100e-9;               % Delay spread: 100ns
Fs = SCS * Nfft;           % SCS = 30 kHz, Nfft = 4096
disp(Fs);
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
rx_grid = zeros(n_Subcarriers, NSlotSymb);  % kích thước giống grid
% Tính độ dài mỗi symbol (bao gồm cả CP)
symbol_lengths = zeros(1, NSlotSymb);
for l = 0:NSlotSymb-1
    if l == 0
        symbol_lengths(l+1) = Nfft + Ncp0;
    else
        symbol_lengths(l+1) = Nfft + Ncp;
    end
end
% Chỉ số bắt đầu mỗi symbol
start_idx = cumsum([1, symbol_lengths(1:end-1)]);
% Demodulation
for l = 0:NSlotSymb-1
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

%% Bước 7: Giải điều chế softbit
rx_data_grid = rx_grid(:, posOFDMforData);       % [3276 × 13]
% Chuẩn bị vector đầu vào cho demodulation
rx_symbols = reshape(rx_data_grid, [], 1);  % [3276*13 x 1] vector

% Tính noise variance (N0) từ SNR
snr_linear = 10^(SNR_dB / 10);
signal_power = mean(abs(rx_symbols).^2);
N0 = signal_power / snr_linear;


% Giải điều chế softbit
softBits = nrSoftModuDemapper(rx_symbols, moduType, noiseVar, 'approx');
rx_bits = softBits < 0;                             % Chuyển sang hardbit
%Dung ham co san
%rx_bits = qamdemod(rx_symbols, 16, 'OutputType', 'bit', 'UnitAveragePower', true);
%rx_bits = reshape(rx_bits, [], 1);  % Đảm bảo giống shape ban đầu của bits_tx
%softBits = qamdemod(rx_symbols, 16, 'OutputType', 'llr', 'UnitAveragePower', true, 'NoiseVariance', noiseVar);
%rx_bits = softBits < 0;
%% Bước 8: Tính MSE và BER
% Bit gốc đã tạo ban đầu
bits_tx = bits;

% So sánh với bit thu
bit_errors = sum(rx_bits(:) ~= bits_tx(:));
total_bits = length(bits_tx);
ber = bit_errors / total_bits;

% Tính MSE giữa symbol gốc và symbol sau giải OFDM
tx_symbols_grid = reshape(modulated_Syms, n_Subcarriers, nbOFDMforData);
mse = mean(abs(rx_data_grid(:) - tx_symbols_grid(:)).^2);

fprintf("BER (bit error rate) = %.6f\n", ber);
fprintf("MSE (mean squared error) = %.6f\n", mse);

figure;
plot(real(tx_signal));
xlabel('Sample Index');
ylabel('Amplitude');
title('Real Part of tx\_signal - Uplink Slot (14 OFDM symbols)');
grid on;

% FFT của tín hiệu phát
Nfft_display = 2^16;  % Zero-padding để có phổ mịn
TX_FFT = fftshift(fft(tx_signal, Nfft_display));
f_axis = linspace(-0.5, 0.5, Nfft_display);  % Tần số chuẩn hóa

figure;
plot(f_axis, 20*log10(abs(TX_FFT)/max(abs(TX_FFT))));
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
title('Spectrum of tx\_signal (14 OFDM symbols)');
grid on;

figure;
imagesc(abs(pusch_grid));   % Lấy trị tuyệt đối vì dữ liệu là số phức
colormap jet;               % Màu sắc dạng phổ
colorbar;                   % Thêm thanh chỉ cường độ màu
xlabel('OFDM Symbol Index (1–14)');
ylabel('Subcarrier Index (1–3276)');
title('Heatmap of |pusch\_grid|');
hold on;
xline(4, '--w', 'LineWidth', 2); % Vạch trắng tại OFDM symbol 4

