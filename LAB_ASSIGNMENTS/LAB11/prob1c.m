clear;
mw=68/60;
cw=4.186;
co=1.9;
tw1=35;
tw2=75;
to1=75;
to2=110;
mo=(mw*cw*(tw2-tw1))/(co*(to2-to1));
dt1=(to1-tw1);
dt2=(to2-tw2);
tlmtd=(dt2-dt1)/log(dt2/dt1);

r=(to2-to1)/(tw2-tw1);
s=(tw2-tw1)/(to2-tw1);
ft=(sqrt(r^2+1)*log((1-s)/(1-r*s)))/((r-1)*log((2-s*(r+1-sqrt(r^2+1)))/(2-s*(r+1+sqrt(r^2+1)))));
U=320;
A=(mw*cw*(tw2-tw1))/(U*ft*tlmtd);
disp("area="+A+"sqm");