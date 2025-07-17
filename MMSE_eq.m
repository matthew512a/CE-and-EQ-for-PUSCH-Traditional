function x_hat = MMSE_eq(y, h, snr_lin)
% MMSE_eq - MMSE Equalization
% Inputs:
%   y: ma trận nhận [n_Subcarriers x num_symbols] (complex)
%   h: vector đáp ứng kênh ước lượng [n_Subcarriers x 1] (complex)
%   snr_lin: SNR theo thang tuyến tính (không phải dB)
% Output:
%   x_hat: ma trận đã cân bằng kênh (complex)

    H_conj = conj(h);
    H_abs2 = abs(h).^2;
    denom = H_abs2 + 1/snr_lin;

    x_hat = (H_conj ./ denom) .* y;
end
