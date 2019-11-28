function [phi_est]=ratchford_rice(T,P,A,B,C)
F=[10 30 5 20 20];%feeed stream
z=F./sum(F);
for i = 1:5
    Psat(i) = antoine_eqn(A(i),B(i),C(i),T);
    k(i) = Psat(i)/P;
end
itr=0;
x=0.1;
f1=0.0;
f0=0.0;
f=0.0;
fprime=0.0;
for i=1:5
    f0=f0+z(i)*k(i);
    f1=f1+z(i)/k(i);
end
if f1>1 & f0 > 1 % checking for function value as f(0) > 0 & f(1) < 0
%newton's method to find phi
    while(1)
        itr=itr+1;
        for i=1:5
            f=f+(z(i)*(k(i)-1))/(1+x*(k(i)-1));
            fprime=fprime-((z(i)*((k(i)-1)^2))/((1+x*(k(i)-1))^2));
        end
        x_1=x-f/fprime;%newton's method equation 
        error=abs(x_1-x); 
        if error<=0.0001
           break;
         elseif itr >= 100
            break;
        end
        x=x_1;
    end
    phi_est=x_1;
elseif f0<1
    phi_est=0;
elseif f1<1
    phi_est=1;
end
end
