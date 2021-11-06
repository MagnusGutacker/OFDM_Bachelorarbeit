function output = remove_cp(data,cp_size,N_sub)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y = length(data)/(N_sub + cp_size) * N_sub;
y = int32(y);
output = zeros(1,y);
for i = 1:length(output)/N_sub
    x1 = i * N_sub - N_sub + 1;
    y1 = i * N_sub;
    x2 = i * (N_sub + cp_size) - N_sub - cp_size + 1;
    y2 = i * (N_sub + cp_size);
    x1 = int32(x1);
    x2 = int32(x2);
    y1 = int32(y1);
    y2 = int32(y2);
%     output(x2 : x2 + cp_size -1) = data(y1 - cp_size + 1 : y1);
%     output(x2 + cp_size : y2) = data (x1 : y1);
    output(x1 : y1) = data (x2 + cp_size: y2);
end
end

