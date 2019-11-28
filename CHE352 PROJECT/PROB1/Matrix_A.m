function [A]=Matrix_A(X,j)
% Matrix - 'A' for jth plate 
Xnew1 = X;
Xnew2 = X;
c = 5;
Xnew1(6,j-1) = X(6,j-1)+(X(6,j-1)/(10^3));% Temperature = Temperature + Delta Temperature
Xnew2(6,j-1) = X(6,j-1)-(X(6,j-1)/(10^3));% Temperature = Temperature - Delta Temperature

%alotting A matrix

A1 = zeros(1,c); %dH/dv (s)
A2 = (Enthalpy_stage(Xnew1,j)-Enthalpy_stage(Xnew2,j))/(2*X(6,j-1)/(10^3));% dH_by_dT = (H + Delta_H)-(H - Delta_H) /2*Delta Temperature
A3 = zeros(1,c);%dH/dl (s)
A4 = zeros(c,c);%dM/dv (s)
A5 = zeros(c,1);%dM/dT (s)
A6 = -eye(c);%dM/dl (s)
A7 = zeros(c,c);%dE/dv (s)
A8 = zeros(c,1);%dE/dT (s)
A9 = zeros(c,c);%dE/dl (s)

for i = (c+2):(2*c+1)
    Xnew1 = X;
    Xnew2 = X;
    Xnew1(i,j-1) = X(i,j-1)+X(i,j-1)/(10^3);
    Xnew2(i,j-1) = X(i,j-1)-X(i,j-1)/(10^3);
    A3(i-6) = (Enthalpy_stage(Xnew1,j)-Enthalpy_stage(Xnew2,j))/(2*X(i,j-1)/(10^3));
end
A =[A1 A2 A3
    A4 A5 A6
    A7 A8 A9];
end