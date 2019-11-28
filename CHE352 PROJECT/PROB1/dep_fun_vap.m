%departure function for vapour phase enthalpies
function enthalpy_departed = dep_fun_vap(T,y)
%critical data taken from NIST database
Pcr = [81 48 47.5 48.9 53.2868]*10^5; %in Pa
Tcr = [513 508 510 562 537]; %in Kelvin
%parameters of Redlich Kwong Equation of State
a = zeros(5,1);
b = zeros(5,1); 
A = zeros(5,1);
B = zeros(5,1);

R = 8.314; %gas constant

for i = 1:5
    a(i) = (0.4278 * R * R * (Tcr(i)^2.5))/Pcr(i);
    A(i) = (a(i)^0.5)/(R*(T^1.25));
    b(i) = 0.0867 * R * Tcr(i) / Pcr(i);
    B(i) = b(i)/(R*T);
end
P = 101325;

AV = 0;
BV = 0;
for i = 1:5
   AV = AV + y(i)*A(i);
   BV = BV + y(i)*B(i);
end
%calculating V through cubic equation
tmp1 = -1*R*T;
tmp2 = AV/(T^0.5) - R*T*BV - P*(BV^2);
tmp3 = -1*(AV*BV)/(T^0.5);
equ = [P tmp1 tmp2 tmp3];
r = roots(equ);
r = r(imag(r)==0);
V = max(r);
z = P*V/R*T;

enthalpy_departed = R*T*(z - 1 - (3*(AV^2)*log(1 + (BV*P)/z)/(2*BV))); %in J/mol
end