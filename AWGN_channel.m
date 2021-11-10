function output = AWGN_channel(data,var)
%AWGN_CHANNEL Summary of this function goes here
%   Detailed explanation goes here

output = awgn(data,var,'measured');

end

