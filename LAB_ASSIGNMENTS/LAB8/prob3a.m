k1=0.004;
k2=k1/4;
y=zeros(3,1);
dydt=@(t,y)[-k1*y(1); k1*y(1)-k2*y(2); k2*y(2)];
[t,y]=ode45(dydt,[0,1200],[100;0;0]);

figure
hold on;
plot(t,y(:,1),'r');
plot(t,y(:,2),'b');
plot(t,y(:,3),'g');
xlabel("time(s)");
ylabel("Concentration of A,B and C (mol/cum)");
legend("concentration of A","concentration of B","concentration of C");
hold off;
figure
hold on;
plot(y(:,2)/62.9933,(100-y(:,1))/100,'r');
xlabel("overall yield of B");
ylabel("conversion of A");
hold off;