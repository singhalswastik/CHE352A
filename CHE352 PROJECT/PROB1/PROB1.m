%Group-8 
% components:- 1-methanol, 2-acetone, 3-methyl acetate, 4-benzene, 5-chloroform
Feed = zeros(5,19);% kmol/h
Feed(:,8) = [10 30 5 20 20];%feed rate
n = 19;% total number of trays with nth being the lowest one
c = 5;% number of components in feed
low_bp = 329.15;%lowest boiling point among all component (acetone)
high_bp = 353.25;%highest boiling point  among all component (benzene)
%reference temperature for enthalpy calculation is taken to be 330 K(feed temperature)
latent_h = Enthalpy_Vapour(Feed(:,8)./85,330);% latent heat of vaporistion for feed
for i = 1:n
    X1(c+1,i) = (low_bp*(n-i)+high_bp*(i-1))/(n-1);% interpolating T for guess
end
%antoine parameters for each component 1-5
A_Antoine = [8.08097 7.11714 7.06524 6.87987 6.95465];
B_Antoine = [1582.271 1210.595 1157.630 1196.760 1170.966];
C_Antoine = [239.726 229.664 219.726 219.161 226.232];

% Rachford Rice to Guess k values
P = 1.013;% in bar
tfeed = (low_bp*(n-8)+high_bp*(8-1))/(n-1); %taking temperature of the feed plate as the flash temperature
phi_guess=ratchford_rice(tfeed,P,A_Antoine,B_Antoine,C_Antoine);
V = phi_guess*85;
L = 85-V;
z=[10 30 5 20 20]./85;
for i = 1:5
    Psat(i) = antoine_eqn(A_Antoine(i),B_Antoine(i),C_Antoine(i),tfeed);
    k(i) = Psat(i)/P;
end
x = zeros(1,c);
y = zeros(1,c);
for i = 1:c
    x(i) = z(i)/((k(i)-1)*phi_guess+1);
    y(i) = k(i)*x(i);
end
l = x.*L;
v = y.*V;
for i = 1:5
    X1(i,:) = v(i);
    X1(i+6,:) = l(i);
end
% Concatenating vapour flow rates of each component, Temperature of each
% stage and liquid flow rates of each component in a single matrix X
X = horzcat(X1(:,1)',X1(:,2)');
for i = 3:19
    X = horzcat(X, X1(:,i)');
end
Tau=1;
Eps=2;
%  Solving Napthali Sandholm algorithm
itr = 1;
while(Tau<Eps)
 H = Enthalpy_Equation(X1);%Enthalpy balance equations at each stage
 M = Material_Equation(X1,Feed);% Material balance equations at each stage
 E = Equilibrium_Equation(X1,A_Antoine,B_Antoine,C_Antoine);% Equilibrium equations at each stage
 
 F = zeros(1,1);%function vector 
 Fprime = zeros((2*c+1),n);% function vector in matrix form
 for j = 1:n
    Fprime(1,j) = H(j);
    for i = 1:c
    Fprime(1+i,j) = M(i,j);
    Fprime(1+c+i,j) = E(i,j);
    end
    F = horzcat(F,Fprime(:,j)');%making one horizontal matrix with columns of fpart
 end
 F = F(2:210);%removing first zero element
 size = 2*c+1;
 dF_by_dx = zeros((size)*n,(size)*n);%derivative of F
 for j = 1:n
    if j == 1
        B = Matrix_B(X1,j);
        C = Matrix_C(X1,j);
        concat = horzcat(B,C);
        dF_by_dx(1+size*(j-1):size*j,1:2*(size)) = concat;
    elseif j == n
        A = Matrix_A(X1,j);
        B = Matrix_B(X1,j);
        concat = horzcat(A,B);
        dF_by_dx(1+size*(j-1):size*j,(size)*n+1-2*(size):(size)*n) = concat;
    else
        A = Matrix_A(X1,j);
        B = Matrix_B(X1,j);
        C = Matrix_C(X1,j);
        concat = horzcat(A,B,C);
        dF_by_dx(1+size*(j-1):size*j,1+size*(j-2):size*(j+1)) = concat;
    end
 end
Eps = (sum(H.^2)+sum(sum(M.^2))+sum(sum(E.^2)))*n*(2*c+1)*(10^-10);
t = Tau;
Tau = (sum(H.^2)/(latent_h^2))+sum(sum(M.^2))+sum(sum(E.^2));
  Xnew = X'-inv(dF_by_dx)*F';% updating X
  X = Xnew';
  fprintf('Current Iteration = %d\n',itr);%number of iterations
 
  itr = itr+1;
  if abs(t-Tau)<0.0001
      break
  end
  for i = 1:n
      X1(:,i) = X(1+(2*c+1)*(i-1):(2*c+1)*i)';
  end
end

%graphs
figure(1);
grid on;
hold on
X1(6,:) = X1(6,:)+(n+1);
plot([1:19],X1(6,:),'b-*');
title('Temperature Profile');
xlabel('Tray Number')
ylabel('Temperature(Kelvin)');

Ls = zeros(1,19);
Vs = zeros(1,19);
for i = 1:19
    Ls(i) = sum(X1(c+2:2*c+1,i));
end
for i = 1:19
    Vs(i) = sum(X1(1:c,i));
end
figure(2)
hold on
grid on
h1 = plot([1:19],X1(1,:)./Vs,'b-*');
h2 = plot([1:19],X1(2,:)./Vs,'g-*');
h3 = plot([1:19],X1(3,:)./Vs,'-*','color',[0.6350, 0.0780, 0.1840]);
h4 = plot([1:19],X1(4,:)./Vs,'-*','color',[0.5 0 0.5]);
h5 = plot([1:19],X1(5,:)./Vs,'-*','color',[0.75 0.75 0.75]);
legend([h1 h2 h3 h4 h5],{'Methanol','Acetone','Methyl Acetate','Benzene','Choloroform'},'Location','southeastoutside');
title('Composition Profile for Vapor Phase');
xlabel('Tray Number')
ylabel('Vapor Molar Flow (kg/h');
figure(3)
hold on
grid on
h1 = plot([1:19],X1(1+c+1,:)./Ls,'b-*');
h2 = plot([1:19],X1(2+c+1,:)./Ls,'g-*');
h3 = plot([1:19],X1(3+c+1,:)./Ls,'-*','color',[0.6350, 0.0780, 0.1840]);
h4 = plot([1:19],X1(4+c+1,:)./Ls,'-*','color',[0.5 0 0.5]);
h5 = plot([1:19],X1(5+c+1,:)./Ls,'-*','color',[0.75 0.75 0.75]);
legend([h1 h2 h3 h4 h5],{'Methanol','Acetone','Methyl Acetate','Benzene','Choloroform'},'Location','southeastoutside');
title('Composition Profile for Liquid Phase');
xlabel('Tray Number')
ylabel('Liquid Molar Flow (kg/h');
hold off