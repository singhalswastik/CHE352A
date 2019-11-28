function [E]=Equilibrium_stage(X,j)
%Equilibrium  equation for the jth tray
%antoine parameters
a = [8.08097 7.11714 7.06524 6.87987 6.95465];
b = [1582.271 1210.595 1157.630 1196.760 1170.966];
C = [239.726 229.664 219.726 219.161 226.232];

E = zeros(5,1);

xv = X(1:5,:); 
xl = X(7:11,:); 
P = 1.013; % pressure in bar

flowv = sum(xv);% total vapor flow rate at each stage
flowl = sum(xl);% total liquid flow rate at each stage

x = (xl(:,j)./flowl(j)); % liquid mole fraction
y = (xv(:,j)./flowv(j)); % vapor mole fraction
gamma=activity(x,X(6,j));% activity coefficients
for i = 1:5
    Psat = antoine_eqn(a(i),b(i),C(i),X(6,j)); 
    E(i)=-xv(i,j)+(Psat/(P*Fugacity_coefficient(X(6,j),i,y)))*gamma(i)*xl(i,j)*flowv(j)/flowl(j);
end
end