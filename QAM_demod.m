function demod = QAM_demod(QAM_signal,size,norm)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

QAM_signal = QAM_signal/norm;

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

