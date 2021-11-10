clc;
clear;
%OFDM
%% Declare parameters

N_sub = 64;                   %Anzahl der Unterträger
symbol_size = 4;                %Symbolgröße z.B. 2 für 4-Qam oder 4 für 16-QAM
signal_length = 100;            %Signallänge in OFDM Symbolen
cp_size = ceil(N_sub/8);        %cyclic prefix Länge                
f = 5e9;                        %Trägerfrequenz
fs = 2e7;                       %Bandbreite
Ts = 1/fs;                      %Symboldauer
var = 20;                       %Varianz/SNR
taps = [0.2,0.1,0.02,0.05,0.05];%Gewichtung der einzelnen Verzögerungen
delays = [1,2,3,4,5];           %Verzögerungen des Kanals

%% generate signal

input_signal = randi([0 1],1,N_sub * signal_length * symbol_size);

%% serial to parallel
%konvertiert den seriellen Datenstrom in N_sub parallele Datenströme 

parallel = serial_to_parallel(input_signal,N_sub,symbol_size);

%% QAM
%Bits werden in QAM-Symbole moduliert

[QAM_modulated,norm]  = QAM(parallel);

%% IFFT
%Die IFFT der parallelen Symbole wird berechnet

for j = 1 : length(QAM_modulated(1,:))
    x = j * length(QAM_modulated(:,1)) - length(QAM_modulated(:,1)) + 1;
    y = j * length(QAM_modulated(:,1));
    ifft_array (x:y) = ifft(QAM_modulated(:,j)); 
end

%% cyclic prefix

cp = cyclic_prefix(ifft_array,cp_size,N_sub);

%% shift to passband

%cp = upconversion(cp,fs,f);

%% Channel
channel_array = cp;
%Tapped Dealy channel 
channel_array = tapped_delay_channel(channel_array,taps,delays);
%AWGN-Channel
channel_array = AWGN_channel(channel_array,var);

%% shift to baseband

%% remove cyclic prefix

no_cp = remove_cp (channel_array,cp_size,N_sub);

%% FFT
%FFT wird von einem OFDM Symbol gebildet, das die Länge von N_sub hat

ifft_array_parallel = serial_to_parallel(no_cp,N_sub,1);

for j = 1 : (signal_length)
    x = j * N_sub - N_sub + 1;
    y = j * N_sub;
    fft_array(x:y) = fft(ifft_array_parallel(:,j)); 
end

%% demodulate QAM
%QAM-Symbole werden in Bits demoduliert

QAM_demodulated = QAM_demod(fft_array,symbol_size,norm);

%% parallel to serial

output_signal = parallel_to_serial(QAM_demodulated);

%% test plots

figure('Name' , "Plots");
hold on;
subplot(2,1,1);
BER = output_signal - input_signal;
plot(1:length(BER),BER);

subplot(2,1,2);
z = interp(ifft_array(1,1:256),4);
plot(-length(z)/2:length(z)/2-1,fftshift(abs(fft(z))));
hold off;


scatterplot(reshape(QAM_modulated,1,[]));

scatterplot(fft_array);


numberOfZeros = sum(BER(:)==0);
1 - numberOfZeros/length(BER)