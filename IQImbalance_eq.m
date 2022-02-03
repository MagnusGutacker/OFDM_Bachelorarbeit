function output = IQImbalance_eq(data,pilot,norm)
%IQ-IMBALANCE_EQ Summary of this function goes here
%   Detailed explanation goes here

data1 = data(1:length(pilot));
data2 = data1-real(pilot);
data3 = exp(-1i * pi/8)*data;
a_est = 1/mean(abs(data3));

output = data3/a_est;
end

