function [ M ] =Material_stage(X,j)
% material eq for the jth tray
c=5;n=19;reflux=9.5;
M=zeros(5,1);
if j==1
    for i=1:5
     M(i)=-1*X(i,2)+X(c+1+i,1)*(1+1/reflux)+X(i,1);
    end
elseif j==19
    for i=1:5
      M(i)=-1*X(c+1+i,18)+X(c+1+i,19)+X(i,19);
    end
else
   for i=1:5 
    M(i)=-1*X(c+1+i,j-1)-1*X(i,j+1)+X(c+1+i,j)+X(i,j);
   end
end
end
