clc;
clear;
%% Declare parameters

signal_in_passband = false;     %true for frequency modulation

N_sub = 64;                     %Anzahl der Unterträger
symbol_size = [2];                %Symbolgröße z.B. 2 für 4-Qam oder 4 für 16-QAM
signal_length = 100;           %Anzahl der OFDM Symbole
cp_size = ceil(N_sub/8);        %cyclic prefix Länge                
fc = 5e9;                       %Trägerfrequenz
fs = 4e7;                       %Bandbreite
SNR = 0:2:16;                 %Varianz/SNR
taps = [1,0.4,0.3];             %Gewichtung der einzelnen Verzögerungen
delays = [0,1,2];               %Verzögerungen des Kanals
ph_off = pi/8;                  %Phasen offset
freq_off = 0.001;                %Frequenz offset (bezogen auf T)
t_off = 0.08;                   %Abtastungs offset
a_off = 1.5;                    %Amplituden offset
IQ_ph_off = pi/8;               %Phasen offset bei IQ_imbalance

for z = 1:length(symbol_size)
pilot_length = 100;
pilot = serial_to_parallel(randi([0 1],1,N_sub * symbol_size(z) * pilot_length),N_sub,symbol_size(z));
pilot = QAM(pilot);

%% generate signal

input_signal = randi([0 1],1,N_sub * signal_length * symbol_size(z));

%% serial to parallel

parallel = serial_to_parallel(input_signal,N_sub,symbol_size(z));

%% QAM

[QAM_modulated,norm]  = QAM(parallel);

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

%% shift to passband

if (signal_in_passband)
    cp = upconversion(cp,fs,fc);
end

pilot_modulated = cp(1:(N_sub + cp_size) * pilot_length );

%% Channel
for i = 1:length(SNR)   %Calculate for each SNR
    %AWGN-Channel
    channel_array = AWGN_channel(cp,SNR(i));

    for a = 1:5         %Calculate for AWGN(a=1) TD(a=2) TD-zf(a=3) TD-MMSE(a=4)
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

        ifft_array_parallel = serial_to_parallel(no_cp,N_sub,1);
        for j = 1 : (signal_length+pilot_length)
            x = j * N_sub - N_sub + 1;
            y = j * N_sub;
            fft_array(x:y) = fft(ifft_array_parallel(:,j));
        end

        fft_array1 = fft_array * (mean(abs(reshape(pilot,1,[])))/mean(abs(fft_array)));
        if (a==5)
            fft_array1 = MLSE(fft_array1,pilot,symbol_size(z),taps,signal_pilot);
        end
        fft_array1 = reshape(fft_array1,1,[]);
        %% Sender-Empfänger Unzulänglichkeiten
        %Phasen Offset
        fft_array_ph_off = phase_offset(fft_array1,ph_off);
        %Frequenz Offset
        fft_array_freq_off = frequency_offset(fft_array_ph_off,freq_off);
        %falsche Abtastung
        %fft_array_off = falsche_abtastung(fft_array_off,t_off);
        %scatterplot(fft_array_off);
        %IQ_imbalance
        fft_array_iq_off = IQ_imbalance(fft_array_freq_off,a_off,IQ_ph_off);

        %% Unzulänglichkeiten Equalizer
        fft_array_freq_eq = frequency_offset_eq(fft_array_freq_off,reshape(pilot,1,[]),freq_off);
        fft_array_iq_eq = IQImbalance_eq(fft_array_freq_eq,reshape(pilot,1,[]),norm);
        fft_array_ph_eq = phase_offset_eq(fft_array_freq_eq,reshape(pilot,1,[]));
        
        %% remove pilot

        no_pilot = fft_array_ph_eq(pilot_length*N_sub+1:end);

        %% demodulate QAM

        QAM_demodulated = QAM_demod(no_pilot,symbol_size(z),norm);

        %% parallel to serial

        output_signal = parallel_to_serial(QAM_demodulated);

        %% Fehlerraten
        %BER
        BitFehler = output_signal - input_signal;
        numberOfZeros = sum(BitFehler(:)==0);
        BER(i,a,z) = 1 - numberOfZeros/length(BitFehler);

        %SER
        SymbolFehler = reshape(QAM_modulated,1,[]) - reshape(QAM(serial_to_parallel(output_signal,N_sub,symbol_size(z))),1,[]);
        numberOfZeros = sum(SymbolFehler(:)==0);
        SER = 1 - numberOfZeros/length(SymbolFehler);
    end
end        
end
BER = reshape(BER,i,a*length(symbol_size));
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
legend({'AWGN-Channel','TD without zf','TD with zf','TD with MMSE','MLSE','AWGN-Channel','TD without zf','TD with zf','TD with MMSE'},'Location','southwest')
% subplot(5,1,5);
% axis([ -1 delays(end)+1 0 1.2]);
% xlabel("Delays");
% ylabel("Taps");
% stem(taps);
%scatterplot(reshape(QAM_modulated,1,[]));
%scatterplot(fft_array);