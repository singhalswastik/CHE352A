function [ M ] = Material_Equation(X,f)
c=5;
n=19;%number of stage
rr=9.5;%reflux ratio
M=zeros(c,n);
% stages in distillation column
for j=2:18
    for i=1:5
    M(i,j)=-1*X(c+1+i,j-1)-1*X(i,j+1)-f(i,j)+X(c+1+i,j)+X(i,j);
    end
end 
% stages as condenser and reboiler 
for i=1:5
    M(i,1)=-1*X(i,2)+X(c+1+i,1)*(1+1/rr)+X(i,1);
    M(i,19)=-1*X(c+1+i,18)+X(c+1+i,19)+X(i,19);
end
end