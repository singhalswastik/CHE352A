function [t]=tsat(a,b,c,p)
    t=(b/(a-log10(p*7.50062)))-c+273;
end