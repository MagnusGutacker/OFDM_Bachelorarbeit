function output = tapped_delay_channel(signal, taps, delays)
%TAPPED_DELY_CHANNEL Summary of this function goes here
%   Detailed explanation goes here

output = zeros(length(signal),1);

for i = 1 : length(signal)-length(delays)
    output(i) = signal(i);
    for j = 1 : length(delays)
        output(i) = output(i) + signal(i+delays(j)) * taps(j);
    end
end
end

