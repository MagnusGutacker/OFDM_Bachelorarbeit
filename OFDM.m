clc;
clear;
%% Declare parameters

signal_in_passband = false;     %true for frequency modulation

N_sub = 64;                     %Anzahl der Unterträger
symbol_size = 2;                %Symbolgröße z.B. 2 für 4-Qam oder 4 für 16-QAM
signal_length = 1000;           %Anzahl der OFDM Symbole
cp_size = ceil(N_sub/8);        %cyclic prefix Länge                
fc = 5e9;                       %Trägerfrequenz
fs = 4e7;                       %Bandbreite
SNR = -10:3:20;                 %Varianz/SNR
taps = [1,0.4,0.3];             %Gewichtung der einzelnen Verzögerungen
delays = [0,1,2];               %Verzögerungen des Kanals

%% generate signal

input_signal = randi([0 1],1,N_sub * signal_length * symbol_size);

%% serial to parallel

parallel = serial_to_parallel(input_signal,N_sub,symbol_size);

%% QAM

[QAM_modulated,norm]  = QAM(parallel);

%% IFFT

for j = 1 : length(QAM_modulated(1,:))
    x = j * length(QAM_modulated(:,1)) - length(QAM_modulated(:,1)) + 1;
    y = j * length(QAM_modulated(:,1));
    ifft_array (x:y) = ifft(QAM_modulated(:,j));
end

%% cyclic prefix

cp = cyclic_prefix(ifft_array,cp_size,N_sub);

%% shift to passband

if (signal_in_passband)
    cp = upconversion(cp,fs,fc);
end
%% Channel
for i = 1:length(SNR)   %Calculate for each SNR
    %AWGN-Channel
    channel_array = AWGN_channel(cp,SNR(i));

    for a = 1:4         %Calculate for AWGN(a=1) TD(a=2) TD-zf(a=3) TD-MMSE(a=4)
        channel_array_out = channel_array;
        if(a==2||a==3||a==4) 
            if (signal_in_passband)
                delays_channel = delays * 50;   %Bei Signal im Passband wird das Signal um den Faktor 50 upgesampled, deshalb müssen die Delays auch 50 mal größer sein
            else
                delays_channel = delays;
            end
            %Tapped-delay-channel
            channel_array_out = tapped_delay_channel(channel_array,taps,delays_channel);
        end
        %% shift to baseband

        baseband_signal = reshape(channel_array_out,[],1);
        if (signal_in_passband)
            baseband_signal = downconversion(baseband_signal,fs,fc);
        end

        %% Equalization 

        if (a==3)           %ZF-Equalization
            equalized = zf_equalizer(baseband_signal,taps,delays);
        elseif (a==4)       %MMSE-Equalization
            equalized = MMSE(baseband_signal,taps,delays,SNR(i));
        else                %No Equalization
            equalized = baseband_signal;
        end

        %% remove cyclic prefix

        no_cp = remove_cp (equalized,cp_size,N_sub);

        %% FFT

        ifft_array_parallel = serial_to_parallel(no_cp,N_sub,1);
        for j = 1 : (signal_length)
            x = j * N_sub - N_sub + 1;
            y = j * N_sub;
            fft_array(x:y) = fft(ifft_array_parallel(:,j));
        end

        %% demodulate QAM

        QAM_demodulated = QAM_demod(fft_array,symbol_size,norm);

        %% parallel to serial

        output_signal = parallel_to_serial(QAM_demodulated);

        %% Fehlerraten
        %BER
        BitFehler = output_signal - input_signal;
        numberOfZeros = sum(BitFehler(:)==0);
        BER(i,a) = 1 - numberOfZeros/length(BitFehler);

        %SER
        SymbolFehler = reshape(QAM_modulated,1,[]) - reshape(QAM(serial_to_parallel(output_signal,N_sub,symbol_size)),1,[]);
        numberOfZeros = sum(SymbolFehler(:)==0);
        SER = 1 - numberOfZeros/length(SymbolFehler);

    end
end        

%% test plots

% if (signal_in_passband)
%     figure('Name' , "Plots");
%     hold on;
%     title('Passband Signal');
%     xlabel('Frequenz f in GHz');
%     ylabel('Amplitude');
%     z = interp(cp,2);
%     f = -(fs*500)/(2):fs*500/(length(z)):(fs*500)/(2)-1;
%     f = f/1000000000; %GHz
%     plot(f,fftshift(abs(fft(z))));
%     hold off;
% end
figure('Name',"BER of equalized and not equalized Signal");
%subplot(5,1,[1 3]);
hold on;
xlabel("SNR in dB");
ylabel("BER");
set(gca,'YScale','log')
axis([ SNR(1) SNR(end) 1/(signal_length*N_sub) 1]);
plot(SNR,BER);
legend({'AWGN-Channel','TD without zf','TD with zf','TD with MMSE'},'Location','southwest')
% subplot(5,1,5);
% axis([ -1 delays(end)+1 0 1.2]);
% xlabel("Delays");
% ylabel("Taps");
% stem(taps);
%scatterplot(reshape(QAM_modulated,1,[]));
%scatterplot(fft_array);