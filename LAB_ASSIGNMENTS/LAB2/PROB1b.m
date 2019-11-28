m=[349.7 0.000 0.000;
349.0 0.022 0.047;
348.2 0.038 0.091;
345.9 0.080 0.206;
343.7 0.122 0.305;
341.7 0.163 0.403;
340.4 0.203 0.449;
338.9 0.237 0.504;
337.5 0.277 0.552;
336.6 0.293 0.587;
334.3 0.399 0.670;
332.1 0.559 0.748;
331.1 0.681 0.792;
330.6 0.819 0.841;
330.7 0.898 0.881;
332.3 1.000 1.000];
tex=m(:,1);
x=m(:,2);
yex=m(:,3);

y=zeros(1,16);
t=zeros(1,16);
tdash=zeros(1,16);
p=94;
t1=tsat(6.9546,1170.966,226.232,p);
t2=tsat(8.11220,1592.864,226.1842,p);
for i=1:16
    x1=x(i);
    x2=(1-x(i));
    t(i)=x1*t1+x2*t2;
    p1=antoine(6.9546,1170.966,226.232,t(i));
    p2=antoine(8.11220,1592.864,226.1842,t(i));
    
    cg12=(58.68/80.67)*exp(-(-268.7676/(1.98709369025*t(i))));
    cg21=(80.67/58.68)*exp(-(1270.3897/(1.98709369025*t(i))));
    
    g1=exp(-log(x1+cg12*x2)+x2*(cg12/(x1+(cg12*x2))-cg21/(x2+(cg21*x1))));
    g2=exp(-log(x2+cg21*x1)-x1*(cg12/(x1+(cg12*x2))-cg21/(x2+(cg21*x1))));
    
    y1=(x1*g1*p1)/((x1*g1*p1)+(x2*g2*p2));
    y2=(x2*g2*p2)/((x1*g1*p1)+(x2*g2*p2));
    while(abs(1-(y1+y2)) > 0.001)
        if((y1+y2)<1)
            t(i)=t(i)+0.01;
        else
            t(i)=t(i)-0.01;
        end
        cg12=(58.68/80.67)*exp(-(-268.7676/(1.98709369025*t(i))));
        cg21=(80.67/58.68)*exp(-(1270.3897/(1.98709369025*t(i))));
    
        g1=exp(-log(x1+cg12*x2)+x2*(cg12/(x1+(cg12*x2))-cg21/(x2+(cg21*x1))));
        g2=exp(-log(x2+cg21*x1)-x1*(cg12/(x1+(cg12*x2))-cg21/(x2+(cg21*x1))));
    
        y1=x1*g1*p1/((x1*g1*p1)+(x2*g2*p2));
        y2=x2*g2*p2/((x1*g1*p1)+(x2*g2*p2));
    end
    y(i)=y1;
    tdash(i)=abs(t(i)-tex(i));
end


figure
plot(x,t);
hold on;
plot(y,t);
plot(x,tex,'-ro');
plot(yex,tex,'-bx');
xlabel('x,y');
ylabel('Temperature(K)');
legend('t-cal','y-calc','t-exp','y-exp');
hold off;


figure
plot(x,y);
hold on;
plot(x,yex,'-ro');
xlabel('x');
ylabel('y');
legend('y-calc','y-exp');
hold off;

figure
plot(x,tdash);
xlabel('x,y');
ylabel('|Texp-Tcal| (K)');
