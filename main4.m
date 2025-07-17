clc; clear; close all;

%% Thông số
n_Subcarriers = 3276;       % 3276
data_Symbols = 13;          % 13
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
pusch_grid = zeros(n_Subcarriers, 14);
DMRS_Signal = readmatrix("dmrs_signal.xlsx", 'Range', 'A1:A3276');
dmrs_col = 4;
pusch_grid(:, dmrs_col) = DMRS_Signal;
data_cols = setdiff(1:14, dmrs_col);
k = 1;
for col = data_cols
    pusch_grid(:, col) = modulated_symbols(k:k + n_Subcarriers - 1);
    k = k + n_Subcarriers;
end

%% OFDM modulation + chèn CP
tx_signal = [];
n_symbol = data_Symbols + pilot_Symbols;
null_left  = floor((Nfft - n_Subcarriers)/2);
null_right = Nfft - null_left - n_Subcarriers;
for l = 0:n_symbol-1
    freq_symbol = pusch_grid(:, l+1);
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

%% So sánh MSE và BER giữa không cân bằng, LS+ZF, MMSE nâng cao
tx_bits = bits;
tx_data_symbols = reshape(modulated_symbols, n_Subcarriers, data_Symbols);

SNR_dB_range = 0:2:30;
mse_noeq = zeros(size(SNR_dB_range));
ber_noeq = zeros(size(SNR_dB_range));
mse_ls = zeros(size(SNR_dB_range));
ber_ls = zeros(size(SNR_dB_range));
mse_mmse = zeros(size(SNR_dB_range));
ber_mmse = zeros(size(SNR_dB_range));

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
    rx_data_grid_ls = rx_data_grid ./ H_ls;
    rx_symbols_ls = reshape(rx_data_grid_ls, [], 1);
    demod_bits_ls = qam16Demodulate3gpp(rx_symbols_ls);
    rx_bits_ls = reshape(demod_bits_ls, n_Subcarriers, data_Symbols, bps);
    bit_errors_ls = sum(rx_bits_ls(:) ~= tx_bits(:));
    ber_ls(idx) = bit_errors_ls / numel(tx_bits);
    mse_ls(idx) = mean(abs(tx_data_symbols(:) - rx_data_grid_ls(:)).^2);

    % MMSE Estimation nâng cao
    noise_var = 1/(10^(SNR/10));
    %H_mmse = (abs(DMRS_Signal).^2) ./ (abs(DMRS_Signal).^2 + noise_var) .* (rx_dmrs ./ DMRS_Signal);
    epsilon = 1e-10;  % nhỏ để tránh chia 0
    H_mmse = (abs(DMRS_Signal).^2) ./ (abs(DMRS_Signal).^2 + noise_var + epsilon) .* ...
         (rx_dmrs ./ (DMRS_Signal + epsilon));
    
    % MMSE Equalization
    W_mmse = conj(H_mmse) ./ (abs(H_mmse).^2 + noise_var);
    rx_data_grid_mmse = rx_data_grid .* W_mmse;
    rx_symbols_mmse = reshape(rx_data_grid_mmse, [], 1);
    demod_bits_mmse = qam16Demodulate3gpp(rx_symbols_mmse);
    rx_bits_mmse = reshape(demod_bits_mmse, n_Subcarriers, data_Symbols, bps);
    bit_errors_mmse = sum(rx_bits_mmse(:) ~= tx_bits(:));
    ber_mmse(idx) = bit_errors_mmse / numel(tx_bits);
    mse_mmse(idx) = mean(abs(tx_data_symbols(:) - rx_data_grid_mmse(:)).^2);

    fprintf("SNR = %2d dB | BER_noeq = %.4f | BER_LS+ZF = %.4f | BER_MMSE = %.4f\n", ...
        SNR, ber_noeq(idx), ber_ls(idx), ber_mmse(idx));
end
disp('Giá trị MSE theo MMSE ở các mức SNR:');
disp(mse_mmse);
%% Vẽ biểu đồ
figure;
plot(SNR_dB_range, mse_noeq, '-o', SNR_dB_range, mse_ls, '-x', SNR_dB_range, mse_mmse, '-s');
legend('MSE không cân bằng', 'MSE LS+ZF', 'MSE MMSE');
xlabel('SNR (dB)'); ylabel('MSE'); title('So sánh MSE giữa các thuật toán cân bằng');

figure;
plot(SNR_dB_range, ber_noeq, '-o', SNR_dB_range, ber_ls, '-x', SNR_dB_range, ber_mmse, '-s');
legend('BER không cân bằng', 'BER LS+ZF', 'BER MMSE');
xlabel('SNR (dB)'); ylabel('BER'); title('So sánh BER giữa các thuật toán cân bằng');

%% Hiển thị giản đồ chòm sao

figure;
subplot(2,2,1);
plot(real(modulated_symbols), imag(modulated_symbols), 'b.');
title('Trước khi truyền (gốc)');
xlabel('In-Phase'); ylabel('Quadrature'); grid on; axis equal;

subplot(2,2,2);
plot(real(rx_symbols_noeq), imag(rx_symbols_noeq), 'r.');
title('Sau truyền - không cân bằng');
xlabel('In-Phase'); ylabel('Quadrature'); grid on; axis equal;

subplot(2,2,3);
plot(real(rx_symbols_ls), imag(rx_symbols_ls), 'g.');
title('Sau LS + ZF');
xlabel('In-Phase'); ylabel('Quadrature'); grid on; axis equal;

subplot(2,2,4);
plot(real(rx_symbols_mmse), imag(rx_symbols_mmse), 'm.');
title('Sau MMSE');
xlabel('In-Phase'); ylabel('Quadrature'); grid on; axis equal;

