function H_interp = interpolate(H_pilot, pilot_loc, Nfft, method)
% Hàm nội suy tín hiệu kênh từ pilot ra toàn bộ dải tần
% H_pilot    : giá trị kênh tại các vị trí pilot
% pilot_loc  : vị trí của các pilot (chỉ mục)
% Nfft       : tổng số subcarrier
% method     : phương pháp nội suy ('linear' hoặc 'spline')

% Xử lý biên: thêm điểm ở đầu nếu cần
if pilot_loc(1) > 1
    pilot_loc = [1, pilot_loc];
    H_pilot = [H_pilot(1), H_pilot];
end

% Xử lý biên: thêm điểm ở cuối nếu cần
if pilot_loc(end) < Nfft
    pilot_loc = [pilot_loc, Nfft];
    H_pilot = [H_pilot, H_pilot(end)];
end

% Nội suy toàn bộ băng thông
H_interp = interp1(pilot_loc, H_pilot, 1:Nfft, method, 'extrap');

end
