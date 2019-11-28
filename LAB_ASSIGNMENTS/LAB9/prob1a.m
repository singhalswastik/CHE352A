clear;
s1=zeros(1,4);
s2=zeros(1,4);
s3=zeros(1,4);
s4=zeros(1,4);
s5=zeros(1,4);
s6=zeros(1,4);
s7=zeros(1,4);
s7dash=zeros(1,4);

s1(1)=100;
s1(2)=300;
s1(4)=1;

for iter=1:1:1000

s2=s1+s7dash;

s3(1)=(1-0.25)*s2(1);
s3(2)=s2(2)-(3*0.25*s2(1));
s3(3)=s2(3)+(2*0.25)*s2(1);
s3(4)=s2(4);

k=[4.8 70 0.051 0.32];
f=sum(s3);
z=s3./f;
phi=rachford(k,z);

V=phi*f;
L=f-V;
for i=1:1:4
    s4(i)=L*(z(i)/(phi*k(i)+1-phi));
    s5(i)=V*((k(i)*z(i))/(phi*k(i)+1-phi));
end

s7=s5./1.01;
s6=s5-s7;

if(abs(s7(2)-s7dash(2))>0.01)
    s7dash=s7;
else
    break;
end
clc;
end
disp("Convergence in "+iter+" steps");
disp("composition of streams:");
disp("columns are of N2, H2, NH3, CO2 respectively");
disp(s1);
disp(s2);
disp(s3);
disp(s4);
disp(s5);
disp(s6);
disp(s7);