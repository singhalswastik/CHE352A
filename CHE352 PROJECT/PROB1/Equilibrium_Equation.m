function [E] = Equilibrium_Equation(X,a,b,C)
c = 5; 
n = 19;%number of stages
E = zeros(c,n);% Equilibrium Equations of each component at each stage
xv = X(1:5,:);% 'xv' represents component flow rate in vapor phase at each stage 
xl = X(7:11,:);% 'xl' represents component flow rate in liquid phase at each stage
P = 1.013;% in bar 
flowv = sum(xv);% total vapor flow rate at each stage
flowl = sum(xl);% total liquid flow rate at each stage

for j = 1:19
    % for wilson equation 
    x = (xl(:,j)./flowl(j));% 'x' is the of each component mole fraction at each stage
    y = (xv(:,j)./flowv(j));% 'y' is the of each component mole fraction at each stage
    gamma=activity(x,X(6,j));
    for i = 1:5
        Psat = antoine_eqn(a(i),b(i),C(i),X(6,j)); % Psat at given temperature 
        E(i,j) = -xv(i,j)+(Psat/(P*Fugacity_coefficient(X(6,j),i,y)))*gamma(i)*(flowv(j)/flowl(j))*xl(i,j);
    end
end
end