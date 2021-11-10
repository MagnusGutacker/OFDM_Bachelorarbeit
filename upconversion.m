function output = upconversion(data,fs,f)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

t = 0:1/fs: (length(data)-1)/fs;
%e = exp(1i*2*pi*f*t);
r = cos(2*pi*f*t);
i = -sin(2*pi*f*t);
output = e .* data;
end

