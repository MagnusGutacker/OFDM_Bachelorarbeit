function Estimation(taps,N_sub,signal_length,symbol_size,snr_start,snr_end,pilot_length)
%MAIN Summary of this function goes here
% Estimation([1,0.3,0.2],64,1000,2,0,20,50,10)
%   Detailed explanation goes here

SNR = snr_start:snr_end;            
delays = [0,1,2];
cp_size = ceil(N_sub/8);

for z = 1:length(symbol_size)

pilot = randi([0 1],1,N_sub * symbol_size(z) * pilot_length);

parallel = reshape(pilot,symbol_size,N_sub,[]);
pilot = permute(parallel,[2 3 1]);

pilot = QAM(pilot);

%% generate signal

input_signal = randi([0 1],1,N_sub * signal_length * symbol_size(z));

%% serial to parallel

parallel = reshape(input_signal,symbol_size,N_sub,[]);
parallel = permute(parallel,[2 3 1]);

%% QAM

[QAM_modulated]  = QAM(parallel);

%% add pilots

signal_pilot = [pilot,QAM_modulated];

%% IFFT

for j = 1 : length(signal_pilot(1,:))
    x = j * length(signal_pilot(:,1)) - length(signal_pilot(:,1)) + 1;
    y = j * length(signal_pilot(:,1));
    ifft_array (x:y) = ifft(signal_pilot(:,j));
end

%% cyclic prefix

cp = cyclic_prefix(ifft_array,cp_size,N_sub);

pilot_modulated = cp(1:(N_sub + cp_size) * pilot_length );

%% Channel
for i = 1:length(SNR)   %Calculate for each SNR
    %AWGN-Channel
    channel_array = AWGN(cp,SNR(i));

    for a = 1:4         %Calculate for AWGN(a=1) TD(a=2) TD-zf(a=3) TD-MMSE(a=4)
        channel_array_out = channel_array;
        if(a==2||a==3||a==4) 
            delays_channel = delays;
            %Tapped-delay-channel
            channel_array_out = tapped_delay_channel(channel_array,taps,delays_channel);
        end
        %% shift to baseband

        baseband_signal = reshape(channel_array_out,[],1);

        %% Equalization 
        H = LS_estimator(pilot_modulated,baseband_signal,SNR(i));

        if (a==3)           %ZF-Equalization
            equalized = zf_equalizer(baseband_signal,H);
        elseif (a==4)       %MMSE-Equalization
            equalized = MMSE(baseband_signal,H,SNR(i));
        else                %No Equalization
            equalized = baseband_signal;
        end

        %% remove cyclic prefix

        no_cp = remove_cp (equalized,cp_size,N_sub);

        %% FFT

        parallel = reshape(no_cp,1,N_sub,[]);
        ifft_array_parallel = permute(parallel,[2 3 1]);
        for j = 1 : (signal_length+pilot_length)
            x = j * N_sub - N_sub + 1;
            y = j * N_sub;
            fft_array(x:y) = fft(ifft_array_parallel(:,j));
        end

        fft_array1 = fft_array * (mean(abs(reshape(pilot,1,[])))/mean(abs(fft_array)));

        fft_array1 = reshape(fft_array1,1,[]);

        
        %% remove pilot

        no_pilot = fft_array1(pilot_length*N_sub+1:end);

        %% demodulate QAM

        QAM_demodulated = QAM_demod(no_pilot,symbol_size(z));

        %% parallel to serial

        output_signal = reshape(QAM_demodulated,1,[]);

        %% Fehlerraten
        %BER
        BitFehler = output_signal - input_signal;
        numberOfZeros = sum(BitFehler(:)==0);
        BER(i,a,z) = 1 - numberOfZeros/length(BitFehler);

    end
end        
end

%% Plot
figure("Name",'Estimated Taps and real Taps');
hold on;
subplot(2,2,1);
axis([ -1 delays(end)+1 0 1.2]);
xlabel("Delays");
ylabel("Taps");
legend({'real Taps'},'Location','northeast')
stem(taps);
subplot(2,2,3);
xlabel("Delays");
ylabel("Taps");
legend({'estimated Taps'},'Location','northeast')
stem(H);
subplot(2,2,[2 4]);
xlabel("SNR in dB");
ylabel("BER");
axis([ SNR(1) SNR(end) 1/(signal_length*N_sub) 1]);
plot(SNR,BER);
set(gca,'YScale','log');
legend({'AWGN-Channel','TD-Channel','ZF-Equalized','MMSE-Equalized'});
end

function [output] = QAM (parallel)

output(1:length(parallel(:,1,1)),1:length(parallel(1,:,1))) = 0;
if length(parallel(1,1,:)) == 1
    for i = 1:length(parallel(1,:))
        for j = 1:length(parallel(:,1))
            if parallel(j,i) == 0
                output(j,i) = -1;
            else
                output(j,i) = 1;
            end
       end
    end
elseif length(parallel(1,1,:)) == 2 %4QAM
    for j = 1:length(parallel(:,1,1))
        for i = 1:length(parallel(1,:,1))
            if parallel(j,i,1) == 0
                output(j,i) = 1;
            elseif parallel(j,i,1) == 1
                output(j,i) = -1;
            end
            if parallel(j,i,2) == 0
                output(j,i) = output(j,i) + 1i;
            elseif parallel(j,i,2) == 1
                output(j,i) = output(j,i) - 1i;
            end
            
        end
    end
    output = output*0.707106781186548;
elseif length(parallel(1,1,:)) == 4 %16QAM
    for j = 1:length(parallel(:,1,1))
        for i = 1:length(parallel(1,:,1))
            I = [parallel(j,i,1),parallel(j,i,4)];
            Q = [parallel(j,i,2),parallel(j,i,3)];
            if I == [0,0]
                output(j,i) = 1;
            elseif I == [0,1]
                output(j,i) = 3;
            elseif I == [1,1]
                output(j,i) = -3;
            elseif I == [1,0]
                output(j,i) = -1;
            end
            if Q == [0,0]
                output(j,i) = output(j,i) + 1i;
            elseif Q == [0,1]
                output(j,i) = output(j,i) + 3i;
            elseif Q == [1,1]
                output(j,i) = output(j,i) - 3i;
            elseif Q == [1,0]
                output(j,i) = output(j,i) - 1i;
            end
        end
    end
    output = output*0.316227766016838;
elseif length(parallel(1,1,:)) == 6 %64QAM
else
    output = parallel;
end
end

function output = tapped_delay_channel(signal, taps, delays)
TD_array = zeros(delays(end)+1,1);
for i = 1:length(delays)
    TD_array(delays(i)+1) = taps(i);
end
output = conv(TD_array,signal);
output(length(signal)+1:end) = [];
end

function demod = QAM_demod(QAM_signal,size)
output = zeros(length(QAM_signal(:,1)),length(QAM_signal(1,:))*size);
counter = 1;
if size == 1
for i = 1:length(QAM_signal)
    if QAM_signal(i) < 0 
        output(i) = 0;
    else
        output(i) = 1;
    end
end
elseif size == 2
    QAM_signal = QAM_signal/0.707106781186548;
    for j = 1 : length(QAM_signal(1,:))
        for i = 1 : length(QAM_signal(:,1))
            I = real(QAM_signal(i,j));
            Q = imag(QAM_signal(i,j));
            if I > 0
                output(counter) = 0;
                counter = counter + 1;
            else
                output(counter) = 1;
                counter = counter + 1;
            end
            if real(Q) > 0 
                output(counter) = 0;
                counter = counter + 1;
            else
                output(counter) = 1;
                counter = counter + 1;
            end
        end
    end
elseif size == 4
    QAM_signal = QAM_signal/0.316227766016838;
    for j = 1 : length(QAM_signal(1,:))
        for i = 1 : length(QAM_signal(:,1))
            tmp = zeros(1:4);
            I = real(QAM_signal(i,j));
            Q = imag(QAM_signal(i,j));
            if I > 2
                tmp(1) = 0;
                tmp(4) = 1;
            elseif I > 0
                tmp(1) = 0;
                tmp(4) = 0;
            elseif I > -2
                tmp(1) = 1;
                tmp(4) = 0;
            else
                tmp(1) = 1;
                tmp(4) = 1;
            end
            if Q > 2
                tmp(2) = 0;
                tmp(3) = 1;
            elseif Q > 0
                tmp(2) = 0;
                tmp(3) = 0;
            elseif Q > -2
                tmp(2) = 1;
                tmp(3) = 0;
            else
                tmp(2) = 1;
                tmp(3) = 1;
            end

            output(counter*4-3:counter*4) = tmp(1:4);
            counter = counter + 1;
        end
    end
elseif size == 6
end
demod = output;
end

function output = zf_equalizer(data,H)
TD_array = reshape(H,[],1);
H = convmtx(TD_array,length(TD_array));
H(length(H(1,:))+1:end,:) = [];
H_inv = H^-1;
f = H_inv(:,1)';

output = conv(f,data);
output(length(data)+1:end) = [];
end

function output = MMSE(data,H_est,SNR)
TD_array = H_est;
for i = 1:length(TD_array)+2
    H(i:i+length(TD_array)-1,i) = TD_array;
end

H_H = conj(H');
SNR = db2mag(SNR);
F = (H_H * H + (1/(SNR))*eye(length(H(1,:))))^(-1)*H_H;
f = F(:,1);

output = conv(f,data);
output(length(data)+1:end) = [];
end

function H = LS_estimator(pilot,data,snr)
% not working on Octave 
% noise_var = 1/(10^(snr/10));
% pilot = reshape(pilot,[],1);
% H_LS = pinv(data(1:length(pilot)))*pilot;
% W = 1/noise_var * autocorr(data(1:length(pilot)),'NumMA',2);
% H = W * H_LS;


H = pinv(diag(pilot))*data(1:length(pilot));
end

function output = cyclic_prefix(data,cp_size,N_sub)
y = length(data)/N_sub * (N_sub+cp_size);
y = int32(y);
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

function output = remove_cp(data,cp_size,N_sub)
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
    output(x1 : y1) = data (x2 + cp_size: y2);
end
end

function output = AWGN (data,snr)
    output = data + sqrt(1/(10^(snr/10))^2).*randn(1,length(data));
end
