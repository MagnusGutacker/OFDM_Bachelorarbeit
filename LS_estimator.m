function H = LS_estimator(pilot,data,snr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
noise_var = 1/(10^(snr/10));
pilot = reshape(pilot,[],1);
H_LS = pinv(data(1:length(pilot)))*pilot;
%data1 = diag(data(1:length(pilot)));
%H_LS = pinv(data1)*pilot;
W = 1/noise_var * autocorr(data(1:length(pilot)),'NumMA',2);
H = W * H_LS;
end

