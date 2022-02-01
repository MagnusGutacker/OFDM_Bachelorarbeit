function output = phase_offset_eq(data,pilot)
%PHASE_PFFSET_EQ Summary of this function goes here
%   Detailed explanation goes here
data1 = data(1:length(pilot));
angles = angle(data1)-angle(pilot);
for i = 1:length(pilot)
if(angles(i)<-2)
    angles(i) = angles(i) + 2*pi;
end
end
ph_off_est = mean(angles);
output = exp(-1i * ph_off_est)*data;


end

