function [p]=antoine(a,b,c,t)
    p=10^(a-(b/(c+t)));
end