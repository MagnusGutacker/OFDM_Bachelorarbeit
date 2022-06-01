function main(taps,N_sub,signal_length,symbol_size,snr_start,snr_end,filename)

% define formats
landscape = "-S930,350";
portrait  = "-S640,480";

% select format
output_format = portrait;

% Create an invisible figure.
fh = figure(1);
set(fh, "visible", "off");

if taps == 1
    taps = [1,0.5+0.5i,0.2+0.3i];
elseif taps == 2
    taps = [1,0.5,0.2];
elseif taps == 3
    taps = [1];
elseif taps == 4
    taps = [1,-0.5+0.3i,0.7-0.6i];
elseif taps == 5
    taps = [1,0.5,0.4,-0.3,0.2,0.1,0.1];
elseif taps == 6
    t = 1:9;
    taps = [1,exp(-t/3).*randi([-1000,1000],1,9)/1000 .* exp(1i*randi([0,2*3*1000],1,9)/1000)];
end

SNR = snr_start:snr_end;            
delays = [0,1,2];

input_signal = randi([0 1],1,N_sub * signal_length * symbol_size);

parallel = reshape(input_signal,symbol_size,N_sub,[]);
parallel = permute(parallel,[2 3 1]);

[QAM_modulated]  = QAM(parallel);

for j = 1 : length(QAM_modulated(1,:))
    x = j * length(QAM_modulated(:,1)) - length(QAM_modulated(:,1)) + 1;
    y = j * length(QAM_modulated(:,1));
    ifft_array (x:y) = ifft(QAM_modulated(:,j));
end


for i = 1:length(SNR)
    channel_array = AWGN(ifft_array,SNR(i));
    for a = 1:2

        channel_array_out = channel_array;
        if a == 2
        channel_array_out = tapped_delay_channel(channel_array,taps,delays);
        end

        parallel = reshape(channel_array_out,1,N_sub,[]);
        ifft_array_parallel = permute(parallel,[2 3 1]);
        for j = 1 : (signal_length)
            x = j * N_sub - N_sub + 1;
            y = j * N_sub;
            fft_array(x:y) = fft(ifft_array_parallel(:,j));
        end

        QAM_demodulated = QAM_demod(fft_array,symbol_size);

        output_signal = reshape(QAM_demodulated,1,[]);

        %BER
        BitFehler = output_signal - input_signal;
        numberOfZeros = sum(BitFehler(:)==0);
        BER(i,a) = 1 - numberOfZeros/length(BitFehler);

    end

end

hold on;

semilogy(SNR,BER);
xlabel("SNR_{S} in dB",'Interpreter','latex');
ylabel("BER");
axis([ SNR(1) SNR(end) 1/(signal_length*N_sub) 1]);
%set(gca,'YScale','log');
legend({'AWGN-Channel','TD-Channel'});

print(filename, "-dpng", output_format);

end

function [output] = QAM (parallel)
%QAM vonverts data to QAM symbols
%   Detailed explanation goes here


% also ppassible to use qammod() function
refconst = 1;
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
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


output = zeros(length(QAM_signal(:,1)),length(QAM_signal(1,:))*size);
counter = 1;
if size == 1
for i = 1:length(QAM_signal)
    if real(QAM_signal(i)) < 0 
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

function output = AWGN (data,snr)
    Es = sum(abs(data).^2)/length(data);
    SNR = 10^(snr/10);
    n = sqrt(Es/(SNR*2))*(randn(1,length(data))+1i*randn(1,length(data)));
    output = data + n;
end