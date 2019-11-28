clear;
E = 10000;
T = 50+273; 
A = 4775;
R = 8.314;
P = 3.5;
k = A*exp(-E/(R*T));
xa=0.94;
vl = 50/17.132 + 50/9.678 + 100/55.5; %total volume of liquid in m3
Ca=(50*(1-xa)/vl); 
r=-k*Ca*Ca;
fa=50000;
V=xa*fa/(-r);
disp("total volume in litres = "+V);