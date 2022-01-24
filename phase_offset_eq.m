function output = phase_offset_eq(data,pilot)
%PHASE_PFFSET_EQ Summary of this function goes here
%   Detailed explanation goes here
data1 = data(1:length(pilot));
ph_off_est = mean(angle(data1)-angle(pilot));
output = exp(-1i * ph_off_est)*data;


end

