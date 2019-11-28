function [C]=Matrix_C(X,j)
% Matrix - 'C' for jth plate
Xnew1 = X;
Xnew2 = X;
c = 5;
Xnew1(6,j+1) = X(6,j+1)+(X(6,j+1)/(10^3));% Temperature = Temperature + Delta Temperature
Xnew2(6,j+1) = X(6,j+1)-(X(6,j+1)/(10^3));% Temperature = Temperature - Delta Temperature

%alotting B matrix

C1 = zeros(1,c);
C2 = (Enthalpy_stage(Xnew1,j)-Enthalpy_stage(Xnew2,j))/(2*X(6,j+1)/(10^3));
C3 = zeros(1,c);
C4 = -eye(c);
C5 = zeros(c,1);
C6 = zeros(c,c);
C7 = zeros(c,c);
C8 = zeros(c,1);
C9 = zeros(c,c);

for i = 1:c
    Xnew1 = X;
    Xnew2 = X;
    Xnew1(i,j+1) = X(i,j+1)+X(i,j+1)/(10^3);
    Xnew2(i,j+1) = X(i,j+1)-X(i,j+1)/(10^3);
    C1(i) = (Enthalpy_stage(Xnew1,j)-Enthalpy_stage(Xnew2,j))/(2*X(i,j+1)/(10^3));
end

for i = 1:c
    Xnew1 = X;
    Xnew2 = X;
     Xnew1(i,j+1) = X(i,j+1)+X(i,j+1)/(10^3);
    Xnew2(i,j+1) = X(i,j+1)-X(i,j+1)/(10^3);
    C7(:,i) = (Equilibrium_stage(Xnew1,j)-Equilibrium_stage(Xnew2,j))./(2*X(i,j+1)/(10^3));
end

 C=[C1 C2 C3
    C4 C5 C6
    C7 C8 C9];

end






