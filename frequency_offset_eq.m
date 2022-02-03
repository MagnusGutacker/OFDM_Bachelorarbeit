function output= frequency_offset_eq(data,pilot,real_off)
%FREQUENCY_OFFSET_EQ Summary of this function goes here
%   Detailed explanation goes here
a_data = angle(data(1:length(pilot)));
a_pilot = angle(pilot);
ph_off_est = a_data-a_pilot;
freq_off_est = mean(ph_off_est(1:end-1)-ph_off_est(2:end));
k = 1:length(data);

output = data.*exp(-1i * real_off * k);
%plot(1:length(output)*10,fftshift(fft(interp(output,10))));
end

