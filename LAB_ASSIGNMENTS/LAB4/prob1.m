z1=0.375;
z2=1-z1;
p=1;
t=90+273;
a1=5.22841;
b1=1807.332;
c1=0.726;
a2=4.07827;
b2=1343.943;
c2=-53.773;
p1=antoine(a1,b1,c1,t);
p2=antoine(a2,b2,c2,t);
k1=p1/p;
k2=p2/p;

if((z1*k1+z2*k2) > 1 && (z1/k1+z2/k2) > 1)
    syms x
    f=@(x)(z1*(k1-1)/(1+x*(k1-1))+z2*(k2-1)/(1+x*(k2-1)));
    f1=@(x)((-(z1*((k1-1)^2)/(1+x*(k1-1)))-(z2*((k2-1)^2)/(1+x*(k2-1)))));
    %fun='670238570036249/(2251799813685248*((3574605706859995*x)/4503599627370496 + 1)) + 2631125873944127/(9007199254740992*((4209801398310603*x)/9007199254740992 - 1))';
    %x=newtonraphson(f,f1,0.5,0.001,50,0.001);
    %while(abs(f(x))>0.001
    x1=0.5;
    xm=x1-(f(x1)/f1(x1));
    ym=f(xm);
    iter=1;
    while (abs(ym) > 0.001) && iter < 50
        iter=iter+1;
        x1=xm;
        xm=x1-(f(x1)/f1(x1));
        ym=f(xm);
        %error(iter)=abs((xm-x1)/x1)*100;
    end
    disp('there exists a unique root for f(phi) : ');
    disp(xm);
    
    x=-0.5:0.01:1.5;
    y=zeros(1,201);
    for i=1:1:201
        y(i)=f(x(i));
    end
    figure
    plot(x,y);
    
    X1=z1/(xm*k1+1-xm);
    X2=z2/(xm*k2+1-xm);
    Y1=k1*X1;
    Y2=k2*X2;
    disp("the liquid composition leaving the drum is [x1,x2]="+X1+","+X2);
    disp("the gaseous composition leaving the drum is [x1,x2]="+Y1+","+Y2);
else
    disp('there does not exist a unique root for f(phi)');
end
