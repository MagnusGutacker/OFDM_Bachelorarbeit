function output = AWGN_channel(data,snr)
%AWGN_CHANNEL Summary of this function goes here
%   Detailed explanation goes here

output = awgn(data,snr,'measured');

end

