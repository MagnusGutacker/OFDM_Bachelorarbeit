function output = plot_tapped_delay(taps,delays)
%PLOT_TAPPED_DELAY Summary of this function goes here
%   Detailed explanation goes here

figure('Name','Tapped-Delay-Channel','NumberTitle','off');
hold on;
stem(0,1,'x');
for i = 1:length(taps)
    stem(delays(i),taps(i),'x');
end
xlim([-1 (length(taps)+1)]);
output = 1;
end

