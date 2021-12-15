function output = zf_equalizer(data,taps,delays)
%ZF_EQUALIZER Summary of this function goes here
%   Detailed explanation goes here

% M = length(taps)+2;
% 
% H = convmtx(taps,M)';
% H(M,:) = [];
% H(1,:) = [];
% 
% I = zeros(size(pinv(H),2), 1,1);
% I(3)=1;

TD_array = zeros(delays(end)+1,1);
for i = 1:length(delays)
    TD_array(delays(i)+1) = taps(i);
end
H = convmtx(TD_array,length(TD_array));
H(length(H(1,:))+1:end,:) = [];
H_inv = H^-1;
f = H_inv(:,1)';

output = conv(f,data);
output(length(data)+1:end) = [];
%output = filter(f,1,data);

end

