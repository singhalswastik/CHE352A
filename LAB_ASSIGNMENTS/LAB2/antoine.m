function [p]=antoine(a,b,c,t)
    p=0.133322*(10^(a-(b/(c+(t-273)))));
end