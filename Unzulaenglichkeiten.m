function Unzulaenglichkeiten(taps,symbol_size,snr_end,nr,ph_off,freq_off,t_off,a_off,IQ_ph_off,filename)
%MAIN Summary of this function goes here
%Unzulaenglichkeiten([1,0.3,0.2],64,1000,2,40,4,pi/8,0.01,0.1,1.5,pi/8)
%   Detailed explanation goes here

% define formats
landscape = "-S930,350";
portrait  = "-S640,480";

% select format
output_format = portrait;

% Create an invisible figure.
fh = figure(1);
set(fh, "visible", "off");

if taps == 1
    taps = [1];
elseif taps == 2
    taps = [1,0.5+0.5i,0.2+0.3i];
elseif taps == 3
    taps = [1,0.5,0.2];
elseif taps == 4
    taps = [1,-0.5+0.3i,0.7-0.6i];
elseif taps == 5
    taps = [1,0.5,0.4,-0.3,0.2,0.1,0.1];
elseif taps == 6
    t = 1:9;
    taps = [1,exp(-t/3).*randi([-1000,1000],1,9)/1000 .* exp(1i*randi([0,2*3*1000],1,9)/1000)];
end

signal_length = 40;
N_sub = 64;
SNR = snr_end;            
delays = 0:length(taps)-1;

input_signal = randi([0 1],1,N_sub * signal_length * symbol_size);

parallel = reshape(input_signal,symbol_size,N_sub,[]);
parallel = permute(parallel,[2 3 1]);

[QAM_modulated]  = QAM(parallel);

for j = 1 : length(QAM_modulated(1,:))
  x = j * length(QAM_modulated(:,1)) - length(QAM_modulated(:,1)) + 1;
  y = j * length(QAM_modulated(:,1));
  ifft_array (x:y) = ifft(QAM_modulated(:,j));
end

channel_array = AWGN(ifft_array,SNR);
channel_array_out = tapped_delay_channel(channel_array,taps,delays);


parallel = reshape(channel_array_out,1,N_sub,[]);
ifft_array_parallel = permute(parallel,[2 3 1]);
for j = 1 : (signal_length)
  x = j * N_sub - N_sub + 1;
  y = j * N_sub;
  fft_array(x:y) = fft(ifft_array_parallel(:,j));
end

%% Sender-Empfänger Unzulänglichkeiten
if nr == 1
  %Phasen Offset
  fft_array_un = phase_offset(fft_array,ph_off);
  a_off = 1;
elseif nr == 2
  %Frequenz Offset
  fft_array_un = frequency_offset(fft_array,freq_off);
  a_off = 1;
elseif nr == 3
  %falsche Abtastung
  fft_array_un = falsche_abtastung(fft_array,t_off);
  a_off = 1;
elseif nr == 4
  %IQ_imbalance
  fft_array_un = IQ_imbalance(fft_array,a_off,IQ_ph_off);
end

x = a_off+0.5;

hold on;
subplot(1,2,1);
plot(real(fft_array_un),imag(fft_array_un),'linestyle','none','marker','.');
title("Mit Unzulaenglichkeit");
axis([-x x -x x],"square");
ylabel("Im");
xlabel("Re");
subplot(1,2,2);
plot(real(fft_array),imag(fft_array),'linestyle','none','marker','.');
title("Ohne Unzulaenglichkeit");
axis([-x x -x x],"square");
ylabel("Im");
xlabel("Re");
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

function output = phase_offset(data,offset)
%offset in Bogenmaß angeben (pi/8 = 22,5°)
output = exp(1i * offset)*data;
end

function output = frequency_offset(data,offset)
% offset in Bogenmaß angegeben. 
k = 1:length(data);
output = data.*exp(1i * offset * k);
end

function output = falsche_abtastung(data,t_off)
%t_off in Abhängigkeit von Ts angegegeben und nur Werte von 0.00 bis 1.00
data = interpft(data,length(data)*100);
output = data(t_off*100:100:end);
end

function output = IQ_imbalance(data,a_imb,phi_imb)
% Imbalance verändert hier nur den Imaginärteil (Q)

Q = imag(data);
Q = Q * exp(1i * phi_imb) * a_imb;
output = real(data) + 1i*Q;
end

function output = AWGN (data,snr)
    Es = sum(abs(data).^2)/length(data);
    SNR = 10^(snr/10);
    n = sqrt(Es/(SNR*2))*(randn(1,length(data))+1i*randn(1,length(data)));
    output = data + n;
end