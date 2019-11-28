function [H]=Enthalpy_Equation(X)
H=zeros(1,19);% to store enthalpy eq
L=zeros(1,19);% to store liquid flow rate
V=zeros(1,19);% to store vapor flow rate
c=5;
B=62; %bottom flow rate
rr=9.5; %reflux ratio
for i=1:19
    L(i)=sum(X(c+2:2*c+1,i));
end
for i=1:19
    V(i)=sum(X(1:c,i));
end
hL=zeros(1,19);%liquid enthalpy
hV=zeros(1,19);% vapor enthalpy

for i=1:19
    hL(i)=Enthalpy_Liquid((X(c+2:2*c+1,i)./L(i)),X(c+1,i));%Enthalpy_liquid(liquid mole fraction,temp)
end

for i=1:19
    hV(i)=Enthalpy_Vapour((X(1:c,i)./V(i)),X(c+1,i)); %Enthalpy_Vapour(Vapour mole fraction,temp)
end
 
   H(1)=L(1)-V(1)*(rr/(rr+1));%condenser
   H(19)=L(19)-B;% reboiler
   
  for j=2:18
    H(j)=-V(j+1)*hV(j+1)-L(j-1)*hL(j-1)+L(j)*hL(j)+V(j)*hV(j);
  end
end
