%vapour phase fugacity coefficient using Redlich Kwong EOS
function phi = Fugacity_coefficient(T,j,y)
%critical data taken from NIST database
Pcr = [81 48 47.5 48.9 53.2868]*10^5; %in Pa
Tcr = [513 508 510 562 537]; %in Kelvin

a = zeros(5,1);
b = zeros(5,1); %parameters of Redlich Kwong Equation of State
A = zeros(5,1);
B = zeros(5,1);
R = 8.314; %gas constant

for i = 1:5
    a(i) = (0.4278 * R * R * (Tcr(i)^2.5))/Pcr(i);
    A(i) = (a(i)^0.5)/(R * (T^1.25));
    b(i) = 0.0867 * R * Tcr(i) / Pcr(i);
    B(i) = b(i)/(R * T);
end
P = 101325;

AV = 0;
BV = 0;
for i = 1:5
    AV = AV + y(i)*A(i);
    BV = BV + y(i)*B(i);
end
tmp1 = (AV^2)*P - BV*P*(1 + BV*P);
tmp2 = -1*((AV*P)^2)*BV;
equ = [1 -1 tmp1 tmp2];
r = roots(equ);

r = r(imag(r)==0);
z = max(r);

term1 = B(j)*(z-1)/BV;
term2 = log(z - BV*P);
term3 = (AV^2)/BV;
term4 = (2*(A(j)/AV) - B(j)/BV);
term5 = log(1 + (BV*P/z));
phi = exp(term1 - term2 - term3*term4*term5);
end
