function output = MMSE(data,taps,SNR)
%MMSE Summary of this function goes here
%   Detailed explanation goes here

% L = length(taps)-1;
% M = length(taps)+1;
% H = zeros(L+M+1,M+1);
% 
% for i = 1:M+1
%     H(i:i+length(taps)-1,i) = taps; 
% end
% H(L+M:L+M+1,:) = [];
% H(1:2,:) = [];
% H_H = conj(H');
% 
% I = zeros(M+1,1);
% I(ceil((M+1)/2))=1;
% 
% SNR = db2mag(SNR);
% f = inv(H_H*H + (1/SNR)*eye(M+1))*H_H*I;
% 
% output = filter(f,1,data);

M = length(taps)+30;

H = convmtx(taps,M)';
H(M-1:M,:) = [];
H(1:2,:) = [];
I = zeros(size(pinv(H),2), 1,1);
I(1)=1;

H_H = conj(H');

SNR = db2mag(SNR);
f = pinv(H_H*H + (1/SNR)*eye(M))*H_H*I;

output = conv(f,data);
output(length(data)+1:end) = [];

end

