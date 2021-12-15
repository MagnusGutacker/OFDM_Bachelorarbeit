function output = MMSE(data,taps,delays,SNR)
%MMSE Summary of this function goes here
%   Detailed explanation goes here

TD_array = zeros(delays(end)+1,1);
for i = 1:length(delays)
    TD_array(delays(i)+1) = taps(i);
end

for i = 1:length(TD_array)+2
    H(i:i+length(TD_array)-1,i) = TD_array;
end
%H(length(H(1,:))+1:end,:) = [];
H_H = conj(H');
SNR = db2mag(SNR);
F = (H_H * H + (1/(SNR))*eye(length(H(1,:))))^(-1)*H_H;
f = F(:,1);

output = conv(f,data);
output(length(data)+1:end) = [];

end

