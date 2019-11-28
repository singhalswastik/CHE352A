function h = Enthalpy_Vapour( y, temp )
%%330 as ref temperature
% comp 1-methanol, 2-acetone, 3-methyl acetate, 4-benzene, 5-chloroform
%heat capacities in the liquid phase
cp_liq{1}=@(T) (105800-362.23.*T+0.9379.*T.^2);
cp_liq{2}=@(T) (135600-177.*T+0.2837.*T.^2+0.000689.*T.^3+0.000689.*T.^4);
cp_liq{3}=@(T) (61260+270.9.*T);
cp_liq{4}=@(T) (162940-344.94.*T+0.85562.*T.^2);
cp_liq{5}=@(T) (124850-166.34.*T+0.43209.*T.^2);
%heat capacities in the vapor phase
cp_vap{1}=@(T) (0.39252*10^5+0.879.*10^5.*((1.9165.*(10^3)./T)./sinh(1.9165.*(10^3)./T)).^2+0.53654.*10^5.*((896.7./T)./sinh(896.7./T)).^2);
cp_vap{2}=@(T) (0.5704*10^5+1.632.*10^5.*((1.607.*(10^3)./T)./sinh(1.607.*10^3./T)).^2+0.968.*10^5.*((731.5./T)./sinh(731.5./T)).^2);
cp_vap{3}=@(T) (0.555*10^5+1.782.*10^5.*((1.26.*(10^3)./T)./sinh(1.26.*10^3./T)).^2+0.853.*10^5.*((562./T)./sinh(562./T)).^2);
cp_vap{4}=@(T) (0.44767*10^5+2.3085.*10^5.*((1.4792.*(10^3)./T)./sinh(1.4792.*10^3./T)).^2+1.6836.*10^5.*((677.66./T)./sinh(677.66./T)).^2);
cp_vap{5}=@(T) (0.3942*10^5+0.6573.*10^5.*((0.928.*(10^3)./T)./sinh(0.928.*10^3./T)).^2+0.493.*10^5.*((399.6./T)./sinh(399.6./T)).^2);

Hvap=[35.21 29.1 30.32 30.72 29.24]; %kJ/mol
Hvap=10^6.*Hvap; % J/Kmol
bp=[64.7 56 57.1 80.1 61.2]; % degree Celsius
bp=bp+273.15; %K


%enthalpy calculation of vapor at given stage
h=0.0;
for i=1:5
    h=h+y(i)*(integral(cp_liq{i},330, bp(i)) + Hvap(i) + integral(cp_vap{i},bp(i), temp)); %J/Kmol
end
h = h + dep_fun_vap(temp,y)*1000;
end
