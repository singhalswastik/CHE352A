function [B]=Matrix_B(X,j)
% Matrix - 'B' for jth plate
Xnew1 = X;
Xnew2 = X;
c = 5;
Xnew1(6,j) = X(6,j)+(X(6,j)/(10^3)); % Temperature = Temperature + Delta Temperature
Xnew2(6,j) = X(6,j)-(X(6,j)/(10^3)); % Temperature = Temperature - Delta Temperature

%alotting B matrix

B1 = zeros(1,c);
B2 = (Enthalpy_stage(Xnew1,j)-Enthalpy_stage(Xnew2,j))/(2*X(6,j)/(10^3));
B3 = zeros(1,c);
B4 = zeros(c,c);
B5 = zeros(c,1);
B6 = zeros(c,c);
B7 = zeros(c,c);
B8 = zeros(c,1);
B9 = zeros(c,c);

for i = 1:c
    Xnew1 = X;
    Xnew2 = X;
    Xnew1(i,j) = X(i,j)+X(i,j)/(10^3);
    Xnew2(i,j) = X(i,j)-X(i,j)/(10^3);
    B1(i) = (Enthalpy_stage(Xnew1,j)-Enthalpy_stage(Xnew2,j))/(2*X(i,j)/(10^3));
end

for i = (c+2):(2*c)+1
    Xnew1 = X;
    Xnew2 = X;
    Xnew1(i,j) = X(i,j)+X(i,j)/(10^3);
    Xnew2(i,j) = X(i,j)-X(i,j)/(10^3);
    B3(i-6) = (Enthalpy_stage(Xnew1,j)-Enthalpy_stage(Xnew2,j))/(2*X(i,j)/(10^3));
end

for i = 1:c
    Xnew1 = X;
    Xnew2 = X;
     Xnew1(i,j) = X(i,j)+X(i,j)/(10^3);
    Xnew2(i,j) = X(i,j)-X(i,j)/(10^3);
    B4(:,i) = (Material_stage(Xnew1,j)-Material_stage(Xnew2,j))./(2*X(i,j)/(10^3));
end
for i = 1:c
    Xnew1 = X;
    Xnew2 = X;
     Xnew1(i,j) = X(i,j)+X(i,j)/(10^3);
    Xnew2(i,j) = X(i,j)-X(i,j)/(10^3);
    B7(:,i) = (Equilibrium_stage(Xnew1,j)-Equilibrium_stage(Xnew2,j))./(2*X(i,j)/(10^3));
end
for i = (c+2):(2*c)+1
    Xnew1 = X;
    Xnew2 = X;
     Xnew1(i,j) = X(i,j)+X(i,j)/(10^3);
    Xnew2(i,j) = X(i,j)-X(i,j)/(10^3);
    B6(:,i-6) = (Material_stage(Xnew1,j)-Material_stage(Xnew2,j))./(2*X(i,j)/(10^3));
end
for i = (c+2):(2*c)+1
    Xnew1 = X;
    Xnew2 = X;
     Xnew1(i,j) = X(i,j)+X(i,j)/(10^3);
    Xnew2(i,j) = X(i,j)-X(i,j)/(10^3);
 B9(:,(i-6)) = (Equilibrium_stage(Xnew1,j)-Equilibrium_stage(Xnew2,j))./(2*X(i,j)/(10^3));
end
 Xnew1 = X;
    Xnew2 = X;
     Xnew1(6,j) = X(6,j)+X(6,j)/(10^3);
    Xnew2(6,j) = X(6,j)-X(6,j)/(10^3);
    B8(:,1) = (Equilibrium_stage(Xnew1,j)-Equilibrium_stage(Xnew2,j))./(2*X(6,j)/(10^3));
    
B= [B1 B2 B3
    B4 B5 B6
    B7 B8 B9];
end