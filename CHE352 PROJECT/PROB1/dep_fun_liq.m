%departure function for liquid phase enthalpies
function enthalpy_departed = dep_fun_liq(T,x)
%critical data taken from NIST database
Pcr = [81 48 47.5 48.9 53.2868]*10^5; %in Pa
Tcr = [513 508 510 562 537]; %in Kelvin
%molar volume from table 13-4
Vm = [40.73; 74.05; 79.84;89.41;80.67].*10^-6; %in m^3/mol

a = zeros(5,1);
b = zeros(5,1); %parameters of Redlich Kwong Equation of State
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
V = 0;
for i = 1:5
    V = V + x(i)*Vm(i);
end

z = P*V/R*T;

AL = 0.0;
BL = 0.0;
for i = 1:5
    AL = AL + x(i)*A(i);
    BL = BL + x(i)*B(i);
end
enthalpy_departed = R*T*(z - 1 - (3*(AL^2)*log(1 + (BL*P)/z)/(2*BL))); %in j/mol
end