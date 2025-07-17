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
BER_zf = zeros(size(snr_range));
MSE_zf = zeros(size(snr_range));
BER_noeq = zeros(size(snr_range));
MSE_noeq = zeros(size(snr_range));

%% Đọc DMRS Signal từ file excell
DMRS_Signal = readmatrix("dmrs_signal.xlsx", 'Range', 'A1:A3276');

for idx = 1:length(snr_range)
    SNR_dB = snr_range(idx);

    %% Tạo dữ liệu & điều chế
    bits = randi([0 1], n_Subcarriers*nbOFDMforData*bps, 1);
    modulated_Syms = nrModuMapper(bits, moduType);

    %% Grid mapping
    pusch_grid = zeros(n_Subcarriers, 14);
    pusch_grid(:, posOFDMforDmrs) = DMRS_Signal;
    k = 1;
    for col = posOFDMforData
        pusch_grid(:, col) = modulated_Syms(k:k + n_Subcarriers - 1);
        k = k + n_Subcarriers;
    end

    %% OFDM modulation
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

    %% Truyền qua kênh
    channel = getChannelInfo(profile, fd, tm, Fs);
    rx_signal = channel(tx_signal);
    rx_signal_awgn = awgn(rx_signal, SNR_dB, 'measured');

    %% OFDM demodulation
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

    %% Trường hợp 1: Không sử dụng ước lượng và cân bằng kênh
    rx_data_grid_noeq = rx_grid(:, posOFDMforData);
    rx_symbols_noeq = reshape(rx_data_grid_noeq, [], 1);
    N0 = mean(abs(rx_symbols_noeq).^2) / (10^(SNR_dB/10));
    softBits_noeq = nrSoftModuDemapper(rx_symbols_noeq, moduType, N0, 'approx');
    rx_bits_noeq = softBits_noeq < 0;
    bit_errors_noeq = sum(rx_bits_noeq(:) ~= bits(:));
    BER_noeq(idx) = bit_errors_noeq / length(bits);
    MSE_noeq(idx) = mean(abs(rx_symbols_noeq - modulated_Syms).^2);

    %% Trường hợp 2: Áp dụng LS/ZF
    H_ls = rx_grid(:, posOFDMforDmrs) ./ DMRS_Signal;
    known_pos = 1:2:3275;
    known_val = H_ls(known_pos);
    all_pos = 1:3276;
    H_interp_freq = interp1(known_pos, known_val, all_pos, 'spline', 'extrap');
    H_est = repmat(H_interp_freq(:), 1, nbOFDMforData);

    rx_data_grid = rx_grid(:, posOFDMforData);
    eq_data = rx_data_grid ./ H_est;

    rx_symbols = reshape(eq_data, [], 1);
    softBits = nrSoftModuDemapper(rx_symbols, moduType, N0, 'approx');
    rx_bits = softBits < 0;
    bit_errors = sum(rx_bits(:) ~= bits(:));
    BER_zf(idx) = bit_errors / length(bits);
    MSE_zf(idx) = mean(abs(rx_symbols - modulated_Syms).^2);
end

%% Vẽ kết quả so sánh
figure;
subplot(2,1,1);
semilogy(snr_range, BER_zf, '-o', snr_range, BER_noeq, '-x'); grid on;
legend('LS+ZF','No Eq'); xlabel('SNR (dB)'); ylabel('BER'); title('So sánh BER');

subplot(2,1,2);
plot(snr_range, MSE_zf, '-o', snr_range, MSE_noeq, '-x'); grid on;
legend('LS+ZF','No Eq'); xlabel('SNR (dB)'); ylabel('MSE'); title('So sánh MSE');

%% Dữ liệu và pilot được vẽ riêng bằng sơ đồ chòm sao I-Q để quan sát trực quan.
figure;
plot(real(modulated_Syms), imag(modulated_Syms), ...
    'o', 'Color', '#233ce6');
axlim = max(max(abs(modulated_Syms))) + 0.05;
ylim([-axlim, axlim]); xlim([-axlim axlim]); axis square;
xlabel("In-phase"); ylabel("Quadrature");
title("Giản đồ chòm sao I-Q trước khi thực hiện điều chế số.")

figure;
plot(real(rx_symbols_noeq), imag(rx_symbols_noeq), ...
    'o', 'Color', '#233ce6');
axlim = max(max(abs(rx_symbols_noeq))) + 0.05;
ylim([-axlim, axlim]); xlim([-axlim axlim]); axis square;
xlabel("In-phase"); ylabel("Quadrature");
title("Giản đồ chòm sao I-Q Khi không sử dụng CE/EQ.")

figure;
plot(real(rx_symbols), imag(rx_symbols), ...
    'o', 'Color', '#233ce6');
axlim = max(max(abs(rx_symbols))) + 0.05;
ylim([-axlim, axlim]); xlim([-axlim axlim]); axis square;
xlabel("In-phase"); ylabel("Quadrature");
title("Giản đồ chòm sao I-Q khi thực hiện LS/ZF.")
