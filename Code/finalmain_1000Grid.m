clc; clear; close all;

%% Cấu hình hệ thống
n_Subcarriers = 3276;
nbOFDMforData = 13;
nbOFDMforDmrs = 1;
u = 1; SCS = 30e3; K = 4;
Nfft = 2048 * K * 2.^(-u);
Ncp0 = 352; Ncp = 288;
posOFDMforDmrs = 4; posOFDMforData = [1 2 3 5 6 7 8 9 10 11 12 13 14];

profile = 'TDL-B'; fd = 25; tm = 100e-9;
Fs = SCS * Nfft;
snr_range = 0:2:30;
nSamples = 1000 ;

DMRS_Signal = readmatrix("dmrs_signal.xlsx", 'Range', 'A1:A3276');

moduTypes = {'QPSK', '16QAM', '64QAM'};
bps_map = containers.Map({'QPSK','16QAM','64QAM'}, [2, 4, 6]);

for m = 1:length(moduTypes)
    moduType = moduTypes{m};
    bps = bps_map(moduType);

    BER_MMSE = zeros(size(snr_range));
    MSE_MMSE = zeros(size(snr_range));
    BER_LS = zeros(size(snr_range));
    MSE_LS = zeros(size(snr_range));

    fprintf('\n===== Đang chạy điều chế %s =====\n', moduType);

    for idx = 1:length(snr_range)
        SNR_dB = snr_range(idx);
        err_mmse = 0; err_ls = 0;
        mse_mmse = 0; mse_ls = 0;

        for i = 1:nSamples
            %% Sinh dữ liệu
            bits = randi([0 1], n_Subcarriers*nbOFDMforData*bps, 1);
            mod_syms = nrModuMapper(bits, moduType);

            %% Ánh xạ vào lưới
            grid = zeros(n_Subcarriers, 14);
            grid(:, posOFDMforDmrs) = DMRS_Signal;
            k = 1;
            for col = posOFDMforData
                grid(:, col) = mod_syms(k:k+n_Subcarriers-1);
                k = k + n_Subcarriers;
            end

            %% OFDM
            tx = [];
            null_left = floor((Nfft - n_Subcarriers)/2);
            null_right = Nfft - null_left - n_Subcarriers;
            for l = 0:13
                f = grid(:,l+1);
                pad = [zeros(null_left,1); f; zeros(null_right,1)];
                t = ifft(ifftshift(pad), Nfft);
                cp = Ncp0*(l==0) + Ncp*(l~=0);
                tx = [tx; t(end-cp+1:end); t];
            end             

            %% Kênh
            ch = getChannelInfo(profile, fd, tm, Fs);
            rx = ch(tx);
            rx = awgn(rx, SNR_dB, 'measured');

            %% Demod OFDM
            rx_grid = zeros(n_Subcarriers,14);
            lens = [Nfft+Ncp0 repmat(Nfft+Ncp,1,13)];
            start = cumsum([1, lens(1:end-1)]);
            for l = 0:13
                seg = rx(start(l+1):start(l+1)+lens(l+1)-1);
                cp = Ncp0*(l==0) + Ncp*(l~=0);
                f = fftshift(fft(seg(cp+1:end), Nfft));
                rx_grid(:,l+1) = f(null_left+1:null_left+n_Subcarriers);
            end

            %% LS + ZF estimation
            H_ls_raw = rx_grid(:, posOFDMforDmrs) ./ DMRS_Signal;
            known = 1:2:3275;
            interp = interp1(known, H_ls_raw(known), 1:3276, 'linear', 'extrap');
            H_est_ls = repmat(interp(:), 1, nbOFDMforData);
            rx_data = rx_grid(:, posOFDMforData);
            eq_ls = rx_data ./ H_est_ls;
            rx_sym_ls = reshape(eq_ls, [], 1);
            sigp = mean(abs(rx_sym_ls).^2);
            N0 = sigp / 10^(SNR_dB/10);
            softBits_LS = nrSoftModuDemapper(rx_sym_ls, moduType, N0, 'approx');
            rx_bits_LS = softBits_LS < 0;
            err_ls = err_ls + sum(rx_bits_LS(:) ~= bits(:));
            mse_ls = mse_ls + mean(abs(rx_sym_ls - mod_syms).^2);

            %% MMSE estimation
            Y_pilot = rx_grid(:, posOFDMforDmrs);
            Xp = DMRS_Signal;
            pilot_loc = find(Xp ~= 0);
            Xp_used = Xp(pilot_loc);
            Yp_used = Y_pilot(pilot_loc);

            H_LS = Yp_used ./ Xp_used;

            tau_rms = 100e-9;
            Delta_f = SCS;
            j2pi_tau_df = 1j * 2 * pi * tau_rms * Delta_f;
            SNR_linear = 10^(SNR_dB / 10);

            Np = length(pilot_loc);
            K1 = repmat(pilot_loc(:), 1, Np);
            K2 = repmat(pilot_loc(:).', Np, 1);
            Rpp = 1 ./ (1 + j2pi_tau_df * (K1 - K2));
            Rpp = Rpp + (1/SNR_linear) * eye(Np);

            K_all = (1:n_Subcarriers)';
            K3 = repmat(K_all, 1, Np);
            K4 = repmat(pilot_loc(:).', n_Subcarriers, 1);
            Rhp = 1 ./ (1 + j2pi_tau_df * (K3 - K4));

            H_MMSE_full = Rhp * (Rpp \ H_LS);
            H_est_time = repmat(H_MMSE_full, 1, nbOFDMforData);
            rx_data_grid = rx_grid(:, posOFDMforData);
            signal_power = mean(abs(rx_data_grid(:)).^2);
            N0 = signal_power / SNR_linear;

            eq_mmse = rx_data_grid .* conj(H_est_time) ./ (abs(H_est_time).^2 + N0);
            rx_sym = reshape(eq_mmse,[],1);
            softBits = nrSoftModuDemapper(rx_sym, moduType, N0, 'approx');
            rx_bits = softBits < 0;
            err_mmse = err_mmse + sum(rx_bits(:) ~= bits(:));
            mse_mmse = mse_mmse + mean(abs(rx_sym - mod_syms).^2);
        end

        BER_MMSE(idx) = err_mmse / (nSamples * length(bits));
        MSE_MMSE(idx) = mse_mmse / nSamples;
        BER_LS(idx) = err_ls / (nSamples * length(bits));
        MSE_LS(idx) = mse_ls / nSamples;

        fprintf("SNR=%2d dB | BER_LS=%.6f | BER_MMSE=%.6f | MSE_LS=%.5f | MSE_MMSE=%.5f\n", ...
            SNR_dB, BER_LS(idx), BER_MMSE(idx), MSE_LS(idx), MSE_MMSE(idx));
    end

    %% Vẽ kết quả
    figure;
    semilogy(snr_range, BER_LS, '-o', snr_range, BER_MMSE, '-s'); grid on;
    xlabel('SNR (dB)'); ylabel('BER');
    legend('LS + ZF', 'MMSE + MMSE');
    title(['BER vs SNR - ', moduType]);

    figure;
    plot(snr_range, MSE_LS, '-o', snr_range, MSE_MMSE, '-s'); grid on;
    xlabel('SNR (dB)'); ylabel('MSE');
    legend('LS + ZF', 'MMSE + MMSE');
    title(['MSE vs SNR - ', moduType]);
end
