function output = MLSE(data,pilot,symbol_size,taps,input)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
M = 2.^symbol_size; %M incoming and M outgoing per state
L = length(taps)-1;

qam = qammod((0:2^(L)-1),2^symbol_size);
norm = modnorm(qam,'avpow',1);
qam = norm*qam;

a = 1; b = 1;
for i = 1:2^M
    state_table(i,:) = [qam(a),qam(b)];
    b = b+1;
    if mod(i,4) == 0
       a = a+1; 
    end
    if b>2^L
        b = 1;
    end
end
Y = zeros(2^M,4);
for i = 1:length(state_table(:,1))
    Y(i,1) = taps*[qam(1),state_table(i,2),state_table(i,1)]';
    Y(i,2) = taps*[qam(2),state_table(i,2),state_table(i,1)]';
    Y(i,3) = taps*[qam(3),state_table(i,2),state_table(i,1)]';
    Y(i,4) = taps*[qam(4),state_table(i,2),state_table(i,1)]';
end

%% Berechnung Branch Metrics
Bm = zeros(M^(L+1),length(data));
for k = 1:length(data)  %Für alle Datenpunkte berechnen
    counter = 1;
    for z = 1:4
        for i = 1:length(state_table(:,1))
            Bm(counter,k)=abs((data(k)-Y(i,z)))^2;
            counter = counter + 1;
        end
    end
end

%% Path Metrics berechnen
Pm = zeros(2^M,length(data)+1);
Pm(2:end,1) = 1000;    

u = 0:M-1;
Index_default = u * M;
for k = 1:length(data)
    Pm_index = Index_default;
    Bm_index = Index_default+1;
    for i = 0:length(Pm(:,1))-1
        if mod(i,4) == 0
            Pm_index = Pm_index + 1;
        end
        if max(Bm_index)>64
            Bm_index = Bm_index - 64;
        end
        possible_Pm = [Pm(Pm_index(1),k)+Bm(Bm_index(1),k),
                       Pm(Pm_index(2),k)+Bm(Bm_index(2),k),
                       Pm(Pm_index(3),k)+Bm(Bm_index(3),k),
                       Pm(Pm_index(4),k)+Bm(Bm_index(4),k)];
        Pm(i+1,k+1) = min(possible_Pm);
        Bm_index = Bm_index + 16;
    end
end

%% Minimalen Pfad finden
for k = 1:length(Pm(1,:))
    o = length(Pm(1,:))-k+1;
    [a,b(o)] = min(Pm(:,o));
end

%% States vom minimalen Pfad auswerten und so eingangs QAM-Symbole erhalten
output = state_table(b,2);
output(end) = [];
output = reshape(output,1,[]); 
end

