function [h]=Enthalpy_stage(X,j)
%enthalpy eq for the jth tray
B = 62; %Bottom flow rate
rr = 9.5; %reflux ratio
c = 5;
L = zeros(1,3); % to store liquid flow rate
V = zeros(1,3); % to store vapor flow rate
hL = zeros(1,3);
hV = zeros(1,3);
 if j == 1
     h = sum(X(7:11,1))-(rr/(rr+1))*sum(X(1:5,1));
 elseif j == 19
     h = sum(X(7:11,19))-B;
 else
     for i = j-1:j+1
       L(i) = sum(X(c+2:2*c+1,i));
       V(i) = sum(X(1:c,i));
       hL(i) = Enthalpy_Liquid((X(c+2:2*c+1,i)./L(i)),X(c+1,i));
       hV(i) = Enthalpy_Vapour((X(1:c,i)./V(i)),X(c+1,i));
     end
     h = -V(j+1)*hV(j+1)-L(j-1)*hL(j-1)+L(j)*hL(j)+V(j)*hV(j);
 end
end
 