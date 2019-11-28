mw=68/60;
cw=4.186;
co=1.9;
tw1=35;
tw2=75;
to1=110;
to2=75;
mo=(mw*cw*(tw2-tw1))/(co*(to1-to2));

dt1=(to1-tw1);
dt2=(to2-tw2);
tlmtd=(dt2-dt1)/log(dt2/dt1);
U=320;
A=(mw*cw*(tw2-tw1))/(U*tlmtd);
disp("we know that the outlet temp of water and oil is same therefore area would be infinite : A="+A);
