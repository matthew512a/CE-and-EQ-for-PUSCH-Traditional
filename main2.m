clc; clear; close all;

%% Thông số
n_Subcarriers = 3276;
data_Symbols = 13;
pilot_Symbols = 1;
bps = 4;
u = 1;
SCS = 30000;
K = 4;
Nfft = 2048 * K * 2.^(-u);
Ncp0 = 352;
Ncp = 288;

%% Bit ngẫu nhiên
bits = randi([0 1], n_Subcarriers, data_Symbols, bps);
bits_reshaped = reshape(bits, [], bps);

%% Điều chế 16QAM
modulated_symbols = qam16Modulate3gpp(bits_reshaped);

%% Lưới tài nguyên và chèn DMRS
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

%% OFDM modulation
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

%% Kênh truyền
delay_profile = 'TDL-B';
fd = 70;
tm = 100e-9;
Fs = SCS * Nfft;
channel = getChannelInfo(delay_profile, fd, tm, Fs);
rx_signal = channel(tx_signal);

%% Khởi tạo biến
SNR_dB_range = 0:2:30;
mse_with_eq = zeros(size(SNR_dB_range));
bler_with_eq = zeros(size(SNR_dB_range));
mse_no_eq = zeros(size(SNR_dB_range));
bler_no_eq = zeros(size(SNR_dB_range));

symbol_lengths = zeros(1, n_symbol);
for l = 0:n_symbol-1
    symbol_lengths(l+1) = Nfft + (Ncp0 * (l==0) + Ncp * (l~=0));
end
start_idx = cumsum([1, symbol_lengths(1:end-1)]);
tx_bits = bits;
tx_data_symbols = reshape(modulated_symbols, n_Subcarriers, data_Symbols);

for idx = 1:length(SNR_dB_range)
    SNR = SNR_dB_range(idx);
    rx_signal_awgn = awgn(rx_signal, SNR, 'measured');
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

    % Trường hợp 1: Không ước lượng/cân bằng
    rx_data_no_eq = reshape(rx_grid(:, data_cols), [], 1);
    bits_no_eq = qam16Demodulate3gpp(rx_data_no_eq);
    rx_bits_no_eq = reshape(bits_no_eq, n_Subcarriers, data_Symbols, bps);
    bler_no_eq(idx) = any(rx_bits_no_eq(:) ~= tx_bits(:));
    mse_no_eq(idx) = mean(abs(tx_data_symbols(:) - reshape(rx_grid(:, data_cols), [], 1)).^2);

    % Trường hợp 2: Có LS + ZF
    rx_dmrs = rx_grid(:, dmrs_col);
    H_ls = rx_dmrs ./ DMRS_Signal;
    rx_data_grid = rx_grid(:, data_cols);
    for l = 1:size(rx_data_grid, 2)
        rx_data_grid(:, l) = rx_data_grid(:, l) ./ H_ls;
    end
    rx_symbols_zf = reshape(rx_data_grid, [], 1);
    demod_bits_zf = qam16Demodulate3gpp(rx_symbols_zf);
    rx_bits_zf = reshape(demod_bits_zf, n_Subcarriers, data_Symbols, bps);
    bler_with_eq(idx) = any(rx_bits_zf(:) ~= tx_bits(:));
    mse_with_eq(idx) = mean(abs(tx_data_symbols(:) - rx_data_grid(:)).^2);

    fprintf("SNR = %2d dB | MSE(noEQ) = %.5f | BLER(noEQ) = %.3f | MSE(EQ) = %.5f | BLER(EQ) = %.3f\n", ...
        SNR, mse_no_eq(idx), bler_no_eq(idx), mse_with_eq(idx), bler_with_eq(idx));
end

%% Vẽ biểu đồ MSE
figure;
plot(SNR_dB_range, mse_no_eq, '-o', 'DisplayName','MSE - Không EQ'); hold on;
plot(SNR_dB_range, mse_with_eq, '-s', 'DisplayName','MSE - Có LS + ZF');
xlabel('SNR (dB)'); ylabel('MSE'); title('So sánh MSE theo SNR'); legend; grid on;

%% Vẽ biểu đồ BLER
figure;
plot(SNR_dB_range, bler_no_eq, '-o', 'DisplayName','BLER - Không EQ'); hold on;
plot(SNR_dB_range, bler_with_eq, '-s', 'DisplayName','BLER - Có LS + ZF');
xlabel('SNR (dB)'); ylabel('BLER'); title('So sánh BLER theo SNR'); legend; grid on;
