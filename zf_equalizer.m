function output = zf_equalizer(data,taps)
%ZF_EQUALIZER Summary of this function goes here
%   Detailed explanation goes here

%L = length(taps)-1;
M = length(taps)+30;
%H = zeros(L+M+1,M+1);

H = convmtx(taps,M);

I = zeros(size(pinv(H),2), 1,1);
I(1)=1;

f = pinv(H)*I;

output = conv(f,data);
output(length(data)+1:end) = [];
%output = filter(f,1,data);

end

