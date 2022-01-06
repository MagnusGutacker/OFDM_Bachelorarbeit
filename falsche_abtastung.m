function output = falsche_abtastung(data,t_off)
%t_off in Abhängigkeit von Ts angegegeben und nur Werte von 0.00 bis 1.00
data = interp(data,100);
output = data(t_off*100:100:end);
end

