function H_MMSE = MMSE_CE(Y, Xp, pilot_loc, Nfft, Nps, h, SNR_dB)
% Hàm ước lượng kênh truyền bằng MMSE
% Y       : tín hiệu nhận tại các vị trí pilot
% Xp      : tín hiệu pilot đã biết
% pilot_loc : chỉ số vị trí của pilot
% Nfft    : số subcarriers
% Nps     : khoảng cách giữa các pilot
% h       : đáp ứng xung kênh (TDL-B)
% SNR_dB  : tỉ số tín hiệu/nhiễu (dB)

SNR = 10^(SNR_dB/10);
H_LS = Y(pilot_loc) ./ Xp;  % Ước lượng LS tại các pilot

% Tính toán tương quan hàm tự động của kênh trong miền tần số
L = length(h);                      % độ dài đáp ứng xung
tau_rms = sqrt(sum((0:L-1).^2 .* abs(h').^2) / sum(abs(h').^2));  % độ trễ RMS
delta_f = 15e3;                     % subcarrier spacing (Hz)
j2pi_tau_df = 1j * 2 * pi * tau_rms * delta_f;

% Tạo chỉ mục tần số
idx_all = 0:(Nfft-1);
idx_pilot = pilot_loc;

% Ma trận tương quan Rhp và Rpp
Rpp = zeros(length(idx_pilot));
Rhp = zeros(Nfft, length(idx_pilot));

for m = 1:length(idx_pilot)
    for n = 1:length(idx_pilot)
        Rpp(m,n) = sinc((idx_pilot(m)-idx_pilot(n))/Nps);  % dạng lý tưởng
    end
end

for m = 1:Nfft
    for n = 1:length(idx_pilot)
        Rhp(m,n) = sinc((idx_all(m)-idx_pilot(n))/Nps);    % dạng lý tưởng
    end
end

% Thêm noise vào Rpp
Rpp = Rpp + (1/SNR)*eye(length(Rpp));

% Tính MMSE estimator
H_MMSE = Rhp * (Rpp \ H_LS(:));  % tương đương với inv(Rpp)*H_LS

end
