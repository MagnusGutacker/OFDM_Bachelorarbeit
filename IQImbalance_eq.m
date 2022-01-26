function output = IQImbalance_eq(data,pilot,norm)
%IQ-IMBALANCE_EQ Summary of this function goes here
%   Detailed explanation goes here

data1 = data(1:length(pilot));
ph_off_est = mean(angle(data1)-angle(pilot));
data = exp(-1i * ph_off_est)*data;
a_est = mean(abs(data))/norm;
Q = imag(data);
data = real(data) + 1i * Q /a_est;

output = data;
end

