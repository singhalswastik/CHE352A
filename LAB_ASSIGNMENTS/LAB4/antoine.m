function [p]=antoine(a,b,c,t)
    p=0.986923*10^(a-(b/(c+t)));
end