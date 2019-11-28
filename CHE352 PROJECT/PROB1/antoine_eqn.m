function [psat]=antoine_eqn(a,b,c,t)
psat=10^(a-(b/(t-273.15+c)));%antoine equation returning pressure in torr and t is in degree K
psat=psat*0.00133322;%psat converted to bar
end