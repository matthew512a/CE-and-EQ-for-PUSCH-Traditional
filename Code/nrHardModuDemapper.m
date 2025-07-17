function bitsOut = nrHardModuDemapper(symbIn, moduType)
% bitsOut = nrHardModuDemapper(symbIn, moduType)
% Giải điều chế cứng các symbol theo sơ đồ điều chế đã chọn (Gray mapping)
% Input:
%   symbIn  - vector symbol đầu vào (phức)
%   moduType - 'BPSK', 'QPSK', '16QAM', '64QAM', '256QAM'
% Output:
%   bitsOut - vector nhị phân (0/1)

switch lower(moduType)
    case 'bpsk'
        symbIn = real(symbIn) ./ (0.5*sqrt(2)); % bỏ chuẩn hóa
        bitsOut = real(symbIn) < 0;
        
    case 'qpsk'
        symbIn = symbIn * sqrt(2); % bỏ chuẩn hóa
        bitsOut = zeros(length(symbIn)*2,1);
        bitsOut(1:2:end) = real(symbIn) < 0;
        bitsOut(2:2:end) = imag(symbIn) < 0;

    case '16qam'
        symbIn = symbIn * sqrt(10);  % bỏ chuẩn hóa
        bitsOut = zeros(length(symbIn)*4,1);
        re = real(symbIn);
        im = imag(symbIn);

        % I-part
        bitsOut(1:4:end) = re > 0;                     % b1
        bitsOut(3:4:end) = abs(re) < 2;                % b3

        % Q-part
        bitsOut(2:4:end) = im > 0;                     % b2
        bitsOut(4:4:end) = abs(im) < 2;                % b4

    case '64qam'
        symbIn = symbIn * sqrt(42);
        bitsOut = zeros(length(symbIn)*6,1);
        re = real(symbIn);
        im = imag(symbIn);

        % Hàm lượng tử mức: [-7 -5 -3 -1 1 3 5 7]
        bitsOut(1:6:end) = re > 0;
        bitsOut(3:6:end) = abs(re) < 4;
        bitsOut(5:6:end) = (mod(round(abs(re)),4) < 2);

        bitsOut(2:6:end) = im > 0;
        bitsOut(4:6:end) = abs(im) < 4;
        bitsOut(6:6:end) = (mod(round(abs(im)),4) < 2);

    case '256qam'
        symbIn = symbIn * sqrt(170);
        bitsOut = zeros(length(symbIn)*8,1);
        re = real(symbIn);
        im = imag(symbIn);

        % I
        bitsOut(1:8:end) = re > 0;
        bitsOut(3:8:end) = abs(re) < 8;
        bitsOut(5:8:end) = mod(floor(abs(re)),8) < 4;
        bitsOut(7:8:end) = mod(floor(abs(re)),4) < 2;

        % Q
        bitsOut(2:8:end) = im > 0;
        bitsOut(4:8:end) = abs(im) < 8;
        bitsOut(6:8:end) = mod(floor(abs(im)),8) < 4;
        bitsOut(8:8:end) = mod(floor(abs(im)),4) < 2;

    otherwise
        error('Unsupported modulation type.');
end

bitsOut = logical(bitsOut); % chuyển sang kiểu logic (0/1)
end
