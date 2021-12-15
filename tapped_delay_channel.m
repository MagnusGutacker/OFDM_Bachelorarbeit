function output = tapped_delay_channel(signal, taps, delays)
%TAPPED_DELY_CHANNEL Summary of this function goes here
%   Detailed explanation goes here

% Verständliche Darstellung:
% for i = 1:length(delays)
%     signal_delayed(i,delays(i)+1:length(signal)+delays(i)) = signal*taps(i);
% end
% 
% output = zeros(1,length(signal)+delays(end));
% 
% for i = 1:length(delays)
%     output = output + signal_delayed(i,:);
% end

TD_array = zeros(delays(end)+1,1);
for i = 1:length(delays)
    TD_array(delays(i)+1) = taps(i);
end
output = conv(TD_array,signal);
output(length(signal)+1:end) = [];

%Falsche Implementierung:
%output = zeros(length(signal),1);
% 
% for i = 1 : length(signal)-length(delays)
%     output(i) = 0;
%     for j = 1 : length(delays)
%         output(i) = output(i) + signal(i+delays(j)) * taps(j);
%     end
% end
end

