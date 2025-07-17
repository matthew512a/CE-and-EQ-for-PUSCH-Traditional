    
clc; clear; close all;

%% Thông số cấu hình 5G
n_Subcarriers = 3276;
nbOFDMforData = 13;
nbOFDMforDmrs = 1;
bps = 6; moduType = '64QAM';
u = 1; SCS = 30e3; K = 4;
Nfft = 2048 * K * 2.^(-u);
Ncp0 = 352; Ncp = 288;
posOFDMforDmrs = 4; posOFDMforData = [1 2 3 5 6 7 8 9 10 11 12 13 14];

%% Thông số mô phỏng kênh truyền
profile = 'TDL-B'; fd = 25; tm = 100e-9;
Fs = SCS * Nfft;

snr_range = 0:2:30;
BER = zeros(size(snr_range));
MSE = zeros(size(snr_range));
BER_LS = zeros(size(snr_range));
MSE_LS = zeros(size(snr_range));
BER_MMSE = zeros(size(snr_range));
MSE_MMSE = zeros(size(snr_range));

%% Đọc DMRS
DMRS_Signal = readmatrix("dmrs_signal.xlsx", 'Range', 'A1:A3276');

for idx = 1:length(snr_range)
    SNR_dB = snr_range(idx);

    %% B1: Sinh dữ liệu và điều chế
    bits = randi([0 1], n_Subcarriers*nbOFDMforData*bps, 1);
    modulated_Syms = nrModuMapper(bits, moduType);

    %% B2: Ánh xạ vào lưới 3276 × 14
    pusch_grid = zeros(n_Subcarriers, 14);
    pusch_grid(:, posOFDMforDmrs) = DMRS_Signal;
    k = 1;
    for col = posOFDMforData
        pusch_grid(:, col) = modulated_Syms(k:k + n_Subcarriers - 1);
        k = k + n_Subcarriers;
    end

    %% B3: OFDM modulation + CP
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

    %% B4: Mô phỏng kênh
    channel = getChannelInfo(profile, fd, tm, Fs);
    rx_signal = channel(tx_signal);
    rx_signal_awgn = awgn(rx_signal, SNR_dB, 'measured');

    impulse = [1; zeros(49,1)];
    h_full = channel(impulse);
    h = h_full (length(h_full));

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

    %% B6: Ước lượng MMSE
    %% B6: Ước lượng kênh MMSE và cân bằng MMSE (theo chuẩn 3GPP)
    Y_pilot = rx_grid(:, posOFDMforDmrs);       % symbol thứ 4 chứa pilot
    Xp = DMRS_Signal;
    pilot_loc = find(Xp ~= 0);                  % chỉ lấy vị trí pilot thực sự
    Xp_used = Xp(pilot_loc);                    % tín hiệu pilot đã biết
    Yp_used = Y_pilot(pilot_loc);               % tín hiệu nhận tại pilot

    H_LS = Yp_used ./ Xp_used;                  % ước lượng LS tại pilot

    %% Thay thế bằng mô hình tương quan MMSE thực sự
    time_MMSE_start = tic;  % Bắt đầu đo thời gian MMSE
    tau_rms = 100e-9;                           % 100 ns
    Delta_f = SCS;                              % Subcarrier spacing = 30kHz
    j2pi_tau_df = 1j * 2 * pi * tau_rms * Delta_f;
    
    SNR_linear = 10^(SNR_dB / 10);

    Np = length(pilot_loc);
    K1 = repmat(pilot_loc(:), 1, Np);
    K2 = repmat(pilot_loc(:).', Np, 1);
    Rpp = 1 ./ (1 + j2pi_tau_df * (K1 - K2));
    Rpp = Rpp + (1/SNR_linear) * eye(Np);       % thêm nhiễu vào ma trận Rpp

    K_all = (1:n_Subcarriers)';                 % toàn bộ subcarriers
    K3 = repmat(K_all, 1, Np);
    K4 = repmat(pilot_loc(:).', n_Subcarriers, 1);
    Rhp = 1 ./ (1 + j2pi_tau_df * (K3 - K4));    % tương quan toàn bộ với pilot

    %% Ước lượng MMSE tại toàn bộ subcarriers
    H_MMSE_full = Rhp * (Rpp \ H_LS);
    if SNR_dB == 20
        writematrix(H_MMSE_full, 'H_MMSE_SNR20dB.txt');
    end

    %% Cân bằng kênh MMSE
    H_est_time = repmat(H_MMSE_full, 1, nbOFDMforData);
    rx_data_grid = rx_grid(:, posOFDMforData);

    signal_power = mean(abs(rx_data_grid(:)).^2);
    N0 = signal_power / SNR_linear;

    eq_data = rx_data_grid .* conj(H_est_time) ./ (abs(H_est_time).^2 + N0);


    %% B7: Giải điều chế softbit và BER
    rx_symbols = reshape(eq_data, [], 1);
    signal_power = mean(abs(rx_symbols).^2);
    N0 = signal_power / (10^(SNR_dB/10));
    softBits = nrSoftModuDemapper(rx_symbols, moduType, N0, 'approx');
    rx_bits = softBits < 0;
    bit_errors = sum(rx_bits(:) ~= bits(:));
    BER_MMSE(idx) = bit_errors / length(bits);
    tx_symbols_grid = reshape(modulated_Syms, n_Subcarriers, nbOFDMforData);
    MSE_MMSE(idx) = mean(abs(eq_data(:) - tx_symbols_grid(:)).^2);

    time_MMSE(idx) = toc(time_MMSE_start);  % Kết thúc đo thời gian MMSE
    %% Không cân bằng
    rx_data_noeq = rx_grid(:, posOFDMforData);
    rx_symbols_noeq = reshape(rx_data_noeq, [], 1);
    signal_power = mean(abs(rx_symbols_noeq).^2);
    N0 = signal_power / (10^(SNR_dB/10));
    softBits_noeq = nrSoftModuDemapper(rx_symbols_noeq, moduType, N0, 'approx');
    rx_bits_noeq = softBits_noeq < 0;
    BER(idx) = sum(rx_bits_noeq(:) ~= bits(:)) / length(bits);
    MSE(idx) = mean(abs(rx_symbols_noeq - modulated_Syms).^2);

    %% LS + ZF
    time_LS_start = tic;  % Bắt đầu đo thời gian LS
    H_ls = rx_grid(:, posOFDMforDmrs) ./ DMRS_Signal;
    known_pos = 1:2:3275;
    known_val = H_ls(known_pos);
    all_pos = 1:3276;
    H_interp_freq = interp1(known_pos, known_val, all_pos, 'spline', 'extrap');
    H_est_time_LS = repmat(H_interp_freq(:), 1, nbOFDMforData);
    rx_data_grid = rx_grid(:, posOFDMforData);
    eq_data_LS = rx_data_grid ./ H_est_time_LS;

    rx_symbols_LS = reshape(eq_data_LS, [], 1);
    signal_power = mean(abs(rx_symbols_LS).^2);
    N0 = signal_power / (10^(SNR_dB/10));
    softBits_LS = nrSoftModuDemapper(rx_symbols_LS, moduType, N0, 'approx');
    rx_bits_LS = softBits_LS < 0;
    BER_LS(idx) = sum(rx_bits_LS(:) ~= bits(:)) / length(bits);
    MSE_LS(idx) = mean(abs(rx_symbols_LS - modulated_Syms).^2);

    time_LS(idx) = toc(time_LS_start);  % Kết thúc đo thời gian LS

    fprintf("SNR = %2d dB | MSE_noeq = %.5f | BER_noeq = %.6f || MSE_LS/ZF = %.5f | BER_LS/ZF = %.6f || MSE_MMSE/MMSE = %.5f | BER_MMSE/MMSE = %.6f\n", ...
        SNR_dB, MSE(idx), BER(idx), MSE_LS(idx), BER_LS(idx), MSE_MMSE(idx), BER_MMSE(idx));

    %% Hiển thị giản đồ chòm sao (chỉ vẽ tại SNR = 20 dB)
        if SNR_dB == 20
            figure;
    
            subplot(2,2,1);
            plot(real(modulated_Syms), imag(modulated_Syms), 'b.');
            title('Trước khi truyền (gốc)');
            xlabel('In-Phase'); ylabel('Quadrature'); grid on; axis equal;
    
            subplot(2,2,2);
            plot(real(rx_symbols_noeq), imag(rx_symbols_noeq), 'r.');
            title('Sau truyền - không cân bằng');
            xlabel('In-Phase'); ylabel('Quadrature'); grid on; axis equal;
    
            subplot(2,2,3);
            plot(real(rx_symbols_LS), imag(rx_symbols_LS), 'g.');
            title('Sau LS + ZF');
            xlabel('In-Phase'); ylabel('Quadrature'); grid on; axis equal;
    
            subplot(2,2,4);
            plot(real(rx_symbols), imag(rx_symbols), 'm.');
            title('Sau MMSE + MMSE');
            xlabel('In-Phase'); ylabel('Quadrature'); grid on; axis equal;
        end
end

%% Vẽ kết quả
figure;
%subplot(2,1,1);
semilogy(snr_range, BER, '-^', snr_range, BER_LS, '-o', snr_range, BER_MMSE, '-s'); grid on;
title('BER vs SNR');
xlabel('SNR (dB)'); ylabel('Bit Error Rate');
legend('Không cân bằng', 'LS + ZF', 'MMSE + MMSE');
grid on;

figure;
%subplot(2,1,2);
plot(snr_range, MSE, '-^', snr_range, MSE_LS, '-o', snr_range, MSE_MMSE, '-s'); grid on;
title('MSE vs SNR');
xlabel('SNR (dB)'); ylabel('Mean Square Error');
legend('Không cân bằng', 'LS + ZF', 'MMSE + MMSE');
grid on;

avg_time_LS = mean(time_LS);
avg_time_MMSE = mean(time_MMSE);

fprintf("\n>>> Thời gian trung bình (LS + ZF): %.4f giây\n", avg_time_LS);
fprintf(">>> Thời gian trung bình (MMSE + MMSE): %.4f giây\n", avg_time_MMSE);

figure;
plot(snr_range, time_LS*1000, '-o', ...
     snr_range, time_MMSE*1000, '-s'); grid on;
title('Thời gian xử lý theo SNR');
xlabel('SNR (dB)'); ylabel('Thời gian (ms)');
legend('LS + ZF', 'MMSE + MMSE');