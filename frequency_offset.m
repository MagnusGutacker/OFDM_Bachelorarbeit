function output = frequency_offset(data,offset)
% offset in Bogenmaß angegeben. 
k = 1:length(data);
output = data.*exp(1i * offset * k);
end

