function output = zf_equalizer(data,taps)
%ZF_EQUALIZER Summary of this function goes here
%   Detailed explanation goes here

% M = length(taps)+4;
% 
% H = convmtx(taps,M)';
% H(M,:) = [];
% H(1,:) = [];
% 
% I = zeros(size(pinv(H),2), 1,1);
% I(3)=1;

f = pinv(taps);

output = conv(f,data);
output(length(data)+1:end) = [];
%output = filter(f,1,data);

end

