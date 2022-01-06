function output = phase_offset(data,offset)
%offset in Bogenmaß angeben (pi/8 = 22,5°)
output = exp(1i * offset)*data;
end

