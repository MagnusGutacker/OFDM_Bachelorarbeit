function output = cyclic_prefix(data,cp_size,N_sub)
%CYCLIC_PREFIX Summary of this function goes here
%   Detailed explanation goes here
y = length(data)/N_sub * (N_sub+cp_size);
output = zeros(1,y);
for i = 1:length(data)/N_sub
    x1 = i * N_sub - N_sub + 1;
    y1 = i * N_sub;
    x2 = i * (N_sub + cp_size) - N_sub - cp_size + 1;
    y2 = i * (N_sub + cp_size);
    output(x2 : x2 + cp_size -1) = data(y1 - cp_size + 1 : y1);
    output(x2 + cp_size : y2) = data (x1 : y1);
end
end

