clc; clear all; close all;

%% Thông số
n_Subcarriers = 3276;
data_Symbols = 13;
pilot_Symbols = 1;
bps = 4;
u = 1;                      % numerology
SCS = 30000;                % Subcarrier Spacing
K = 4;
Nfft = 2048 * K * 2.^(-u);  % Số mẫu IFFT theo chuẩn 3GPP (4096)
Ncp0 = 352;                 % tại l = 0
Ncp = 288;                  % tại l = 1:13

%% Tạo bit ngẫu nhiên và reshape
bits = randi([0 1], n_Subcarriers, data_Symbols, bps);     % [3276 x 13 x 4]
bits_reshaped = reshape(bits, [], bps);                    % [3276*13 x 4]

%% Điều chế 16QAM - 3GPP
modulated_symbols = qam16Modulate3gpp(bits_reshaped);      % [3276*13 x 1]

%% Tạo lưới tài nguyên và chèn DMRS
grid = zeros(n_Subcarriers, 14);
DMRS_Signal = readmatrix("dmrs_signal.xlsx", 'Range', 'A1:A3276');
dmrs_col = 4;
grid(:, dmrs_col) = DMRS_Signal;
data_cols = setdiff(1:14, dmrs_col);
k = 1;
for col = data_cols
    grid(:, col) = modulated_symbols(k:k + n_Subcarriers - 1);
    k = k + n_Subcarriers;
end

%% OFDM modulation + chèn CP
tx_signal = [];
n_symbol = data_Symbols + pilot_Symbols;
null_left  = floor((Nfft - n_Subcarriers)/2);
null_right = Nfft - null_left - n_Subcarriers;
for l = 0:n_symbol-1
    freq_symbol = grid(:, l+1);
    freq_padded = [zeros(null_left,1); freq_symbol; zeros(null_right,1)];
    time_symbol = ifft(ifftshift(freq_padded), Nfft);
    cp_len = Ncp0 * (l==0) + Ncp * (l~=0);
    cp = time_symbol(end - cp_len + 1:end);
    tx_symbol = [cp; time_symbol];
    tx_signal = [tx_signal; tx_symbol];
end

fprintf("Tín hiệu phát cuối cùng có %d mẫu\n", length(tx_signal));

%% Kênh truyền
delay_profile = 'TDL-B';
fd = 70;
tm = 100e-9;
Fs = SCS * Nfft;
channel = getChannelInfo(delay_profile, fd, tm, Fs);
rx_signal = channel(tx_signal);

%% So sánh MSE và BER giữa có và không dùng cân bằng
tx_bits = bits;
tx_data_symbols = reshape(modulated_symbols, n_Subcarriers, data_Symbols);

SNR_dB_range = 0:2:30;
mse_noeq = zeros(size(SNR_dB_range));
ber_noeq = zeros(size(SNR_dB_range));
mse_eq = zeros(size(SNR_dB_range));
ber_eq = zeros(size(SNR_dB_range));

symbol_lengths = zeros(1, n_symbol);
for l = 0:n_symbol-1
    symbol_lengths(l+1) = Nfft + (Ncp0 * (l==0) + Ncp * (l~=0));
end
start_idx = cumsum([1, symbol_lengths(1:end-1)]);

for idx = 1:length(SNR_dB_range)
    SNR = SNR_dB_range(idx);
    rx_signal_awgn = awgn(rx_signal, SNR, 'measured');

    % OFDM demodulation
    rx_grid = zeros(n_Subcarriers, n_symbol);
    for l = 0:n_symbol-1
        idx_start = start_idx(l+1);
        idx_end = idx_start + symbol_lengths(l+1) - 1;
        rx_symbol_with_cp = rx_signal_awgn(idx_start:idx_end);
        cp_len = Ncp0 * (l==0) + Ncp * (l~=0);
        rx_symbol = rx_symbol_with_cp(cp_len+1:end);
        freq_symbol = fftshift(fft(rx_symbol, Nfft));
        rx_grid(:, l+1) = freq_symbol(null_left+1 : null_left+n_Subcarriers);
    end

    % No Equalization
    rx_data_grid = rx_grid(:, data_cols);
    rx_symbols_noeq = reshape(rx_data_grid, [], 1);
    demod_bits_noeq = qam16Demodulate3gpp(rx_symbols_noeq);
    rx_bits_noeq = reshape(demod_bits_noeq, n_Subcarriers, data_Symbols, bps);
    bit_errors_noeq = sum(rx_bits_noeq(:) ~= tx_bits(:));
    ber_noeq(idx) = bit_errors_noeq / numel(tx_bits);
    mse_noeq(idx) = mean(abs(tx_data_symbols(:) - rx_data_grid(:)).^2);

    % LS Estimation + ZF Equalization
    rx_dmrs = rx_grid(:, dmrs_col);
    H_ls = rx_dmrs ./ DMRS_Signal;
    rx_data_grid_zf = rx_data_grid ./ H_ls;
    rx_symbols_eq = reshape(rx_data_grid_zf, [], 1);
    demod_bits_eq = qam16Demodulate3gpp(rx_symbols_eq);
    rx_bits_eq = reshape(demod_bits_eq, n_Subcarriers, data_Symbols, bps);
    bit_errors_eq = sum(rx_bits_eq(:) ~= tx_bits(:));
    ber_eq(idx) = bit_errors_eq / numel(tx_bits);
    mse_eq(idx) = mean(abs(tx_data_symbols(:) - rx_data_grid_zf(:)).^2);

    fprintf("SNR = %2d dB | MSE_noeq = %.5f | BER_noeq = %.6f || MSE_eq = %.5f | BER_eq = %.6f\n", ...
        SNR, mse_noeq(idx), ber_noeq(idx), mse_eq(idx), ber_eq(idx));
end

%% Vẽ biểu đồ
figure;
plot(SNR_dB_range, mse_noeq, '-o', SNR_dB_range, mse_eq, '-x');
legend('MSE không cân bằng', 'MSE LS+ZF');
xlabel('SNR (dB)'); ylabel('MSE'); title('So sánh MSE giữa có và không cân bằng kênh');

figure;
plot(SNR_dB_range, ber_noeq, '-o', SNR_dB_range, ber_eq, '-x');
legend('BER không cân bằng', 'BER LS+ZF');
xlabel('SNR (dB)'); ylabel('BER'); title('So sánh BER giữa có và không cân bằng kênh');