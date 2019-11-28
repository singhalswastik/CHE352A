clear;
mc=3.8;
mh=1.9;

tc1=38;
tc2=55;
th2=93;
th1=th2-(mc*(tc2-tc1)/mh);
r=(th2-th1)/(tc2-tc1);
s=(tc2-tc1)/(th2-tc1);
dt1=(th1-tc1);
dt2=(th2-tc2);
tlmtd=(dt2-dt1)/log(dt2/dt1);
U=1419;


r=0.95/100;
m=pi*r*r*0.366*1000;  %mass flow rate per tube
nt=ceil(mc/m);       %number of tubes per pass
disp("number of tubes per pass : "+nt);

%case 1 one tube pass
A1=(mc*4186*(tc2-tc1))/(U*tlmtd);
L1=A1/(2*pi*r*nt);

%case 2 two tube pass
ft2=0.89;
np=2;
A2=(mc*4186*(tc2-tc1))/(U*ft2*tlmtd);
L2=A2/(2*pi*r*nt*np);
disp("two tube passes will  be required");
disp("length required : "+L2+"m");
disp("total area of contact : "+A2+"sqm");
disp("the results are same as obtained in aspen");
