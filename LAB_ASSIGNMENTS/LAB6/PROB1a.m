a=9;
b=-6;
y=@(x)((9*x)/(1+8*x))-(0.6)*x*(1-x);
x=0:0.01:1;
yeq=zeros(1,101);
for i=1:1:101
    yeq(i)=y(x(i));
end
xb=0.001;
xd=0.999;

m=@(x)(9/((1+8*x)^2))-(0.6)*(1-2*x);

f=@(x)(m(x)-(y(x)-xd)/(x-xd));
xc=fsolve(f,0.990);

xf=0.4;
yf=m(xc)*(xf-xd)+xd;          %q=1;
mb=(yf-xb)/(xf-xb);         %slope of bottom operating line
mt=m(xc);         %slope of top operating line
r=mt/(1-mt);
s=1/(mb-1);
disp("minimum reflux ratio= "+r);
disp("minimum reboiler ratio= "+s);
%{
figure
plot(x,x);
hold on;
plot(x,yeq);
xlabel('x');
ylabel('y');
plot([xf; xb], [yf; xb]);
plot([xf; xf], [xf; yf]);
plot([xf; xd], [yf; xd]);
hold off;
%}

r=1.5*r;            %r=1.5*rmin
mt=r/(r+1);
yt=@(x)((mt*x)+(xd/(r+1)));        %equation of top operating line
yf=yt(xf);
mb=(yf-xb)/(xf-xb);
s=1/(mb-1);
xbl=@(y)(y+xb/s)/mb;            %equation of bottom operating line
xtl=@(y)(y-(xd/(r+1)))/mt; 
figure
plot(x,x);
hold on;
plot(x,yeq);
xlabel('x');
ylabel('y');
plot([xf; xb], [yf; xb]);
plot([xf; xf], [xf; yf]);
plot([xf; xd], [yf; xd]);

x1=xb;
y1=y(xb);
py1=x1;
c=0;
while x1<xf
    plot([x1; x1], [py1; y1]);
    plot([x1; xbl(y1)], [y1; y1]);
    c=c+1;
    py1=y1;
    x1=xbl(y1);
    y1=y(x1);
end
plot([x1; xtl(py1)], [py1; py1]);
x1=xtl(py1);
y1=y(x1);
while x1<xd
    plot([x1; x1], [py1; y1]);
    plot([x1; xtl(y1)], [y1; y1]);
    c=c+1;
    py1=y1;
    x1=xtl(y1);
    y1=y(x1);
end
hold off;
disp("number of stages with 1.5rmin= "+(c));
disp("Extra two stages for condenser and reboiler will be required");
clear;