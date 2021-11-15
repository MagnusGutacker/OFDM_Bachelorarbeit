function output = upconversion(data,fs,fc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

data = interp(data,50);
t = 0:1/(fs*50):(length(data)-1)/(fs*50);

output = data .* exp(-1i*2*pi*fc*t);

end

