y=@(x)(4.64*x/((3.64)*x+1));    %equilibrium equation
x=0:0.01:1;
yeq=zeros(1,101);
for i=1:1:101
    yeq(i)=y(x(i));
end
figure
plot(x,x);
hold on;
plot(x,yeq);
xlabel('x');
ylabel('y');
xb=0.04;
xd=0.96;
x1=xb;
c=0;
while x1<xd                  %to find the minimum number of stages
    y1=y(x1);
    plot([x1; x1], [x1; y1]);
    plot([x1; y1], [y1; y1]);
    c=c+1;
    x1=y1;
end
hold off;
disp("minimum number of stages= "+(c));

xf=0.45;
yf=y(xf);          %q=1;
mb=(yf-xb)/(xf-xb);         %slope of bottom operating line
mt=(yf-xd)/(xf-xd);         %slope of top operating line
r=mt/(1-mt);
s=1/(mb-1);
disp("minimum reflux ratio= "+r);
disp("minimum reboiler ratio= "+s);

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

r=1.5*r;            %r=1.5*rmin
mt=r/(r+1);
yt=@(x)((mt*x)+(xd/(r+1)));        %equation of top operating line
yf=yt(xf);
mb=(yf-xb)/(xf-xb);
s=1/(mb-1);
xbl=@(y)(y+xb/s)/mb;            %equation of bottom operating line

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
xs=@(y)(y/(4.64-(3.64*y)));
y1=xd;
x1=xs(xd);
px1=x1;
plot([x1; y1], [y1; y1]);
c=c+1;
while x1>xf
    plot([x1; x1], [yt(x1); y1]);
    c=c+1;
    px1=x1;
    y1=yt(x1);
    x1=xs(y1);
    plot([x1; px1], [y1; y1]);
end
hold off;
disp("number of stages with 1.5rmin= "+(c));
disp("Extra two stages for condenser and reboiler will be required");