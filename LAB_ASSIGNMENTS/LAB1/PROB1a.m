m = [4.618 0.13970 0.61830;
5.809 0.20550 0.71900;
7.050 0.27480 0.78790;
8.289 0.34420 0.83620;
9.539 0.41440 0.87240;
10.812 0.48680 0.90080;
12.022 0.55590 0.92240;
13.248 0.62650 0.94040;
14.512 0.69890 0.95580;
14.684 0.70870 0.95770;
15.692 0.76710 0.96810;
16.849 0.83360 0.97850;
16.989 0.84130 0.97970;
17.744 0.88520 0.98580];

pex=m(:,1);
x=m(:,2);
yex=m(:,3);
t=308-273;
p1= 10^(6.87987-(1196.76/(219.161+t)));
p2= 10^(6.99053-(1453.430/(215.310+t)));
y=zeros(1,14);
p=zeros(1,14);
pdash=zeros(1,14);

for i=1:1:14
    y(i)=x(i)*p1/(x(i)*p1+(1-x(i))*p2);
    p(i)=(x(i)*p1*0.133322)/y(i);
    pdash(i)=abs(p(i)-pex(i));
end
figure
plot(x,p);
hold on;
plot(y,p);
plot(x,pex,'-ro');
plot(yex,pex,'-bx');
xlabel('x,y');
ylabel('pressure(kPa)');
legend('x-cal','y-calc','x-exp','y-exp');
hold off;



figure
plot(x,y);
hold on;
plot(x,yex);
xlabel('x');
ylabel('y');
legend('y-calc','y-exp');
hold off;

figure
plot(x,pdash);
xlabel('x,y');
ylabel('|Pcal-Pexp| (kPa)');