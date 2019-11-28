function [h] = Enthalpy_Liquid( x, temp )
% comp 1-methanol, 2-acetone, 3-methyl acetate, 4-benzene, 5-chloroform
%heat capacities in the liquid phase
%ref temp 330
cp_liq{1}=@(T) (105800-362.23.*T+0.9379.*T.^2);
cp_liq{2}=@(T) (135600-177.*T+0.2837.*T.^2+0.000689.*T.^3+0.000689.*T.^4);
cp_liq{3}=@(T) (61260+270.9.*T);
cp_liq{4}=@(T) (162940-344.94.*T+0.85562.*T.^2);
cp_liq{5}=@(T) (124850-166.34.*T+0.43209.*T.^2);
h=0.0;
%enthalpy calculation of liquid at given stage
for i=1:5
    h=h+x(i)*(integral(cp_liq{i},330, temp)); %J/kmol
end
h = h + dep_fun_liq(temp,x)*10^(3);
end