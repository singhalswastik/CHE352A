function [gamma] = activity(x,T)
%%% activity coeff
lambda = [       0          0.6643*10^3    0.8346*10^3    1.6795*10^3    1.7026*10^3;
            -0.2150*10^3         0        -0.0652*10^3    0.4949*10^3   -0.0732*10^3;
            -0.0984*10^3    0.1613*10^3         0         0.0102*10^3    0.0791*10^3;
             0.2161*10^3   -0.1679*10^3    0.2005*10^3         0         0.1484*10^3;
            -0.3723*10^3   -0.3131*10^3   -0.4314*10^3   -0.2085*10^3         0       ];
RU = 1.987;% gas constant value in cal/mol*Kelvin
gamma = zeros(1,5);
lngamma=zeros(1,5);
A = zeros(5,5);
vl = [40.73; 74.05; 79.84;89.41;80.67].*10^-3; %in m^3/mol
for i = 1:5
    for j = 1:5
        A(i,j) = (vl(j)/vl(i))*exp(-lambda(i,j)/(RU*T));
    end
end
first_term = zeros(1,5);
for i = 1:5
 for j = 1:5
     first_term(i) = first_term(i)+x(j,1)*A(i,j);
 end
end
for  i = 1:5
    sec_term = 0.0;
    for j = 1:5
        sec_term = sec_term+(x(j,1)*A(j,i)/first_term(j));
    end
    lngamma(i) = 1-log(first_term(i))-sec_term;
end
gamma = exp(lngamma);
end

