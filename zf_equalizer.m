function output = zf_equalizer(data,taps)
%ZF_EQUALIZER Summary of this function goes here
%   Detailed explanation goes here

H = zeros(length(taps)+length(taps)-1,length(taps));

for i = 1:length(taps)
    H(i:i+length(taps)-1,i) = taps; 
end

I = zeros(length(taps)+length(taps)-1,1);
I(1) = 1;

f = linsolve(H,I);

output = filter(f,1,data);

end

