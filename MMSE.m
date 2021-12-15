function output = MMSE(data,taps,delays,SNR)
%MMSE Summary of this function goes here
%   Detailed explanation goes here

TD_array = zeros(delays(end)+1,1);
for i = 1:length(delays)
    TD_array(delays(i)+1) = taps(i);
end
% H = convmtx(TD_array,length(TD_array));
% H(length(H(1,:))+1:end,:) = [];
for i = 1:length(TD_array)+2
    H(i:i+length(TD_array)-1,i) = TD_array;
%     if(i+length(TD_array)-1>length(TD_array)+2)
%         
%     end
end

H(length(H(1,:))+1:end,:) = [];
H_H = conj(H');
SNR = db2mag(SNR);
F = (H_H*H + (1/(SNR))*eye(length(H(1,:))))^(-1)*H_H;
f = F(1,:);

% M = length(taps)+2;
% 
% H = convmtx(taps,M)';
% H(7,:) = [];
% H(1,:) = [];
% I = zeros(size(pinv(H),2), 1,1);
% I(1)=2;
% 
% H_H = conj(transpose(H));
% 
% SNR = db2mag(SNR);
% f = (H*H_H + (1/SNR)*eye(M))^(-1)*H_H;
% 
% SNR = db2mag(SNR);
% f = taps./(abs(taps).^2+SNR);

output = conv(f,data);
output(length(data)+1:end) = [];

end

