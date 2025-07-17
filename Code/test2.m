clc; clear; close all;

%% Thông số cấu hình 5G
n_Subcarriers = 3276;               % Number of Subcarriers
nbOFDMforData = 13;                 % Number of OFDM Symbol for Data
nbOFDMforDmrs = 1;                  % Number of OFDM Symbol for Pilot
bps = 4; moduType = '16QAM';        % bits per symbol 
u = 1; SCS = 30e3; K = 4;
Nfft = 2048 * K * 2.^(-u);          % Số mẫu IFFT theo chuẩn 3GPP (4096)
Ncp0 = 352; Ncp = 288;              % Number of CP symbols at l = 0, l = 1:13
posOFDMforDmrs = 4; posOFDMforData = [1 2 3 5 6 7 8 9 10 11 12 13 14]; 

%% Thông số mô phỏng kênh truyền
profile = 'TDL-B';                  % Delay Profile
fd = 25;                            % Maximum Doppler Shift
tm = 100e-9;                        % Delay Spread
Fs = SCS * Nfft;                    % Sample Rate

%% Tính toán BER, MSE.
snr_range = 0:2:30;
BER = zeros(size(snr_range));
MSE = zeros(size(snr_range));
BER_LS = zeros(size(snr_range));
MSE_LS = zeros(size(snr_range));
%% Đọc file DMRS Signal
DMRS_Signal = readmatrix("dmrs_signal.xlsx", 'Range', 'A1:A3276');

for idx = 1:length(snr_range)
    SNR_dB = snr_range(idx);
    
    %% B1: Sinh dữ liệu và điều chế
    bits = randi([0 1], n_Subcarriers*nbOFDMforData*bps, 1);        % [3276*13*4 x 1]
    modulated_Syms = nrModuMapper(bits, moduType);

    %% B2: Ánh xạ vào lưới 3276 * 14
    pusch_grid = zeros(n_Subcarriers, 14);
    pusch_grid(:, posOFDMforDmrs) = DMRS_Signal;
    k = 1;
    for col = posOFDMforData
        pusch_grid(:, col) = modulated_Syms(k:k + n_Subcarriers - 1);
        k = k + n_Subcarriers;
    end

    %% B3: OFDM modulation + thêm CP
    tx_signal = [];
    n_symbol = nbOFDMforData + nbOFDMforDmrs;
    null_left  = floor((Nfft - n_Subcarriers)/2);
    null_right = Nfft - null_left - n_Subcarriers;
    for l=0:n_symbol-1
        freq_symbol = pusch_grid(:, l+1);
        freq_padded = [zeros(null_left,1); freq_symbol; zeros(null_right,1)];
        time_symbol = ifft(ifftshift(freq_padded), Nfft);
        cp_len = Ncp0 * (l==0) + Ncp * (l~=0);
        cp = time_symbol(end - cp_len + 1:end);
        tx_symbol = [cp; time_symbol];
        tx_signal = [tx_signal; tx_symbol];
    end

    %% B4: Mô phỏng kênh truyền
    channel = getChannelInfo(profile, fd, tm, Fs);
    rx_signal = channel(tx_signal);
    rx_signal_awgn = awgn(rx_signal, SNR_dB, 'measured');
    
    impulse = [1; zeros(49,1)];
    h_full = channel(impulse);
    h = h_full(1:50);   % hoặc length(h_full) nếu không cố định L

    %% B5: OFDM demodulation
    rx_grid = zeros(n_Subcarriers, n_symbol);
    symbol_lengths = [Nfft+Ncp0 repmat(Nfft+Ncp,1,n_symbol-1)];
    start_idx = cumsum([1, symbol_lengths(1:end-1)]);
    for l = 0:n_symbol-1
        idx_start = start_idx(l+1);
        idx_end = idx_start + symbol_lengths(l+1) - 1;
        rx_symbol_with_cp = rx_signal_awgn(idx_start:idx_end);
        cp_len = Ncp0 * (l==0) + Ncp * (l~=0);
        rx_symbol = rx_symbol_with_cp(cp_len+1:end);
        freq_symbol = fftshift(fft(rx_symbol, Nfft));
        rx_grid(:, l+1) = freq_symbol(null_left+1 : null_left+n_Subcarriers);
    end
    
    %% B6: Ước lượng kênh LS và cân bằng ZF
    % Step 6a: LS estimation tại các subcarrier có DMRS
    H_ls = rx_grid(:, posOFDMforDmrs) ./ DMRS_Signal;
    % Step 6b: Nội suy spline theo tần số (chỉ dùng subcarrier lẻ)
    known_pos = 1:2:3275;
    known_val = H_ls(known_pos);
    all_pos = 1:3276;
    H_interp_freq = interp1(known_pos, known_val, all_pos, 'spline', 'extrap');
    if SNR_dB == 20
        writematrix(H_interp_freq, 'H_LS_SNR20dB.txt');
    end
    % Step 6c: Nội suy theo thời gian (giả định kênh không đổi)
    H_est_time = repmat(H_interp_freq(:), 1, nbOFDMforData);

    % Step 6d: Cân bằng ZF
    rx_data_grid = rx_grid(:, posOFDMforData);
    eq_data = rx_data_grid ./ H_est_time;
    
    % No Equalization (use for compare)
    % No Equalization (dùng trực tiếp rx_data_grid)
    rx_data_noeq = rx_grid(:, posOFDMforData);
    rx_symbols_noeq = reshape(rx_data_noeq, [], 1);
    signal_power = mean(abs(rx_symbols_noeq).^2);
    N0 = signal_power / (10^(SNR_dB/10));
    softBits_noeq = nrSoftModuDemapper(rx_symbols_noeq, moduType, N0, 'approx');
    rx_bits_noeq = softBits_noeq < 0;
    BER(idx) = sum(rx_bits_noeq(:) ~= bits(:)) / length(bits);
    MSE(idx) = mean(abs(rx_symbols_noeq - modulated_Syms).^2);

    %% B7: Giải điều chế softbit và tính BER (with LS/ZF)
    rx_symbols = reshape(eq_data, [], 1);
    signal_power = mean(abs(rx_symbols).^2);
    N0 = signal_power / (10^(SNR_dB/10));
    softBits = nrSoftModuDemapper(rx_symbols, moduType, N0, 'approx');
    rx_bits = softBits < 0;

    bit_errors = sum(rx_bits(:) ~= bits(:));
    BER_LS(idx) = bit_errors / length(bits);

    tx_symbols_grid = reshape(modulated_Syms, n_Subcarriers, nbOFDMforData);
    MSE_LS(idx) = mean(abs(eq_data(:) - tx_symbols_grid(:)).^2);

    fprintf("SNR = %2d dB | MSE_noeq = %.5f | BER_noeq = %.6f || MSE_eq = %.5f | BER_eq = %.6f\n", ...
        SNR_dB, MSE(idx), BER(idx), MSE_LS(idx), BER_LS(idx));
end

%% B8: Vẽ kết quả
figure;
%subplot(2,1,1);
semilogy(snr_range, BER_LS, '-o',snr_range, BER, '-*'); grid on;
title('BER vs SNR '); xlabel('SNR (dB)'); ylabel('Bit Error Rate');
legend('LS + ZF', 'Không ước lượng/cân bằng');
%subplot(2,1,2);
figure;
plot(snr_range, MSE_LS, '-o', snr_range, MSE, '-*'); grid on;
title('MSE vs SNR '); xlabel('SNR (dB)'); ylabel('Mean Square Error');
legend('LS + ZF', 'Không ước lượng/cân bằng');