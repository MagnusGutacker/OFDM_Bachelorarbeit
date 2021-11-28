function output = MMSE(data,taps,SNR)
%MMSE Summary of this function goes here
%   Detailed explanation goes here

L = length(taps)-1;
M = length(taps)+2;
H = zeros(L+M+1,M+1);

for i = 1:M+1
    H(i:i+length(taps)-1,i) = taps; 
end
H(M-1:M+1,:) = [];
H(1:2,:) = [];

H_H = conj(H');

I = zeros(M+1,1);
I(ceil(M+1/2))=1;

SNR = db2mag(SNR);
f = (H_H*H + 1/SNR+eye(M+1))^-1*H_H*I;

output = filter(f,1,data);

end

