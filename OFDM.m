clc;
clear;

%OFDM

%% Declare parameters

N_sub = 256;         %Anzahl der Unterträger
FFT_size = 256;
symbol_size = 4;    %Symbolgröße z.B. 2 für 4-Qam oder 4 für 16-QAM
f = 10000;          %Frequenz
fs = 100;           %Sampling rate
signal_length = 1000;%Signallänge in OFDM Symbolen
cp_size = 25;%cyclic prefix Länge


%% generate signal

input_signal = randi([0 1],1,N_sub * signal_length);

%% serial to parallel
%konvertiert den seriellen Datenstrom in N_sub parallele Datenströme 

parallel = serial_to_parallel(input_signal,N_sub,symbol_size);

%% QAM
%Bits werden in QAM-Symbole moduliert

QAM_modulated  = QAM(parallel);

%QAM_modulated(201:256,:) = 0;

%% IFFT
%Die IFFT der parallelen Symbole wird berechnet (länge der IFFT also N_sub)

for j = 1 : length(QAM_modulated(1,:))
    x = j * length(QAM_modulated(:,1)) - length(QAM_modulated(:,1)) + 1;
    y = j * length(QAM_modulated(:,1));
    ifft_array (x:y) = ifft(QAM_modulated(:,j)); 
end

ifft_array = (ifft_array); %abs?
%% cyclic prefix

cp = cyclic_prefix(ifft_array,cp_size,FFT_size);

%% shift to passband


%% shift to baseband


%% remove cyclic prefix

no_cp = remove_cp (cp,cp_size,FFT_size);

%% FFT
%FFT wird von einem OFDM Symbol gebildet, das die Länge von N_sub hat
ifft_array2 = serial_to_parallel(no_cp,N_sub,1);
fft_array = zeros(FFT_size, signal_length/symbol_size , 'double') + 1i*zeros(FFT_size, signal_length/symbol_size , 'double');
for j = 1 : (signal_length/symbol_size)
    x = j * N_sub - N_sub + 1;
    y = j * N_sub;
    fft_array(x:y) = fft(ifft_array2(:,j)); 
end

%fft_array(201:256,:) = [];

%% demodulate QAM
%QAM-Symbole werden in Bits demoduliert

QAM_demodulated = QAM_demod(fft_array,symbol_size);

%% parallel to serial

output_signal = parallel_to_serial(QAM_demodulated);

%% test plots
t = 0:1/fs:(length(parallel_to_serial(ifft_array))-1)/fs;
f = -length(ifft_array)/2 : 1 : length(ifft_array)/2 - 1;
figure('Name' , "Plots");
hold on;
subplot(2,1,1);
plot(t,ifft_array);
subplot(2,1,2);
plot(-128:127,fftshift(fft(ifft_array(1:256))));
%plot(1:256,QAM_modulated(:,1));
%f = fftshift(fourier_trans);

%plot(t,abs(fftshift(parallel_to_serial(ifft_array))));
%plot(1:50000,fftshift(fourier_trans));
BER = output_signal - input_signal;
%plot(1:64,input_signal(1:64),1:64,output_signal(1:64));
%plot(1:16000,BER);
%plot(1:1000 , parallel);
%plot(1:32000 , real(inverse_fourier_trans));
%figure;
%plot(1:17750 , fft(cp));
%plot(1:64,parallel_data(:,1),1:64,20*real(QAM_modulated(:,1)));