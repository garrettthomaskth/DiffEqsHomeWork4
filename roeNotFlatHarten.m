function [Q,x,t,cons] = roeNotFlatHarten(xSteps, tSteps, fix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Chosen Variable
L=10; %L long
T=15;
%For entropy fix
d = 0.001;

%Given Variables
H=1.05;
%w=0.1*L;
%alpha=H/5;
g=9.8;

r=L/6;
B0 = H/10;

dx=L/xSteps;
dt=T/tSteps;
%Initialize Matrix
Q = zeros(xSteps+2, 2*(tSteps+1));
%Initial Conditions
Q(1:(end),1) = H;
for j = 1:xSteps+2
    x = (j-1.5)*dx;
        if abs(x - L/2) < r            
             Q(j,1) = Q(j,1)-B0*cos(pi*(x-L/2)/(2*r))^2;
        end
end
%H+alpha*exp(-((dx/2:dx:(L-dx/2)) - L/2).^2/w^2);
% i=1;
% Q(1,2*i-1) = Q(2,2*i-1);
% Q(1,2*i) = -Q(2,2*i);
% Q(end,2*i-1) = Q(end-1,2*i-1);
% Q(end,2*i) = -Q(end-1,2*i);
%Choose initial conditions for momentum
Q(:,2) = (Q(:,1)).*sqrt(g*(Q(:,1)))-0.01;
M = Q(1,2);

f = @(u) [ u(2) , u(2)^2./u(1) + 0.5*g*u(1).^2];
 
%Ffun = @(u1,u2) 0.5*(f(u1)+f(u2)) - 0.5*(abs(lambda1)*W1 + abs(lambda2)*W2);

F = zeros(xSteps+1,2);
S = zeros(xSteps+2,2);
for i = 1:tSteps+1
    % Ghost point values
    Q(1,2*i-1) = H;
    Q(1,2*i) = M;
    Q(end,2*i-1:2*i) = 2*Q(end-1,2*i-1:2*i) - Q(end-2,2*i-1:2*i);
    Q(end,2*i-1) = H;
    Q(end,2*i) = M;
    
    for j = 1:xSteps+1
        cHat = sqrt(0.5*g*(Q(j,2*i-1)+Q(j+1,2*i-1)));
        uHat = (sqrt(Q(j,2*i-1))*Q(j,2*i)/Q(j,2*i-1)+sqrt(Q(j+1,2*i-1))*Q(j+1,2*i)/Q(j+1,2*i-1))/(sqrt(Q(j,2*i-1))+sqrt(Q(j+1,2*i-1)));
        
        delta = Q(j+1,2*i-1:2*i)-Q(j,2*i-1:2*i);
        
        W1 = ((uHat+cHat)*delta(1)-delta(2))/(2*cHat) * [1, uHat - cHat];
        W2 = (-(uHat-cHat)*delta(1)+delta(2))/(2*cHat) * [1, uHat + cHat];
                
        lambda1 = uHat - cHat;
        lambda2 = uHat + cHat;
        
        phiLambda1 = lambda1;
        phiLambda2 = lambda2;
        
        % Phi Delta function for the lambdas 
        if fix
            phiLambda1 = phiDelta(lambda1, d);
            phiLambda2 = phiDelta(lambda2, d);
        end
        F(j,:) = 0.5*(f(Q(j+1,(2*i-1):(2*i)))+f(Q(j,(2*i-1):(2*i)))) - 0.5*(abs(phiLambda1)*W1 + abs(phiLambda2)*W2);
        
        x = (j-1.5)*dx;
        
        if abs(x - L/2) < r            
            S(j,2) = g*(Q(j,2*i-1))*pi*B0*sin(pi*(x-L/2)/r)/(2*r);
        end
    end
    for j = 2:xSteps+1
        Q(j,2*i+1:2*i+2) = Q(j,2*i-1:2*i) - dt/dx * (F(j,:)-F(j-1,:))+ dt*S(j,:);
    end
end
M = Q(1:end,2*(1:tSteps+1));

Q(:,2*(1:tSteps+1)-1) = Q(:,2*(1:tSteps+1)-1)+B*ones(1,tSteps+1);

x = linspace(-dx,L+dx,xSteps+2);
t = linspace(0,T,tSteps+1);

Q = Q(1:end,2*(1:tSteps+1)-1);

figure(2)
mesh(t,x,Q)
title('H')
rotate3d on
figure(3)
mesh(t,x,M./(Q-B*ones(1,tSteps+1)))
title('V')
rotate3d on
figure(4)
plot(dx:dx:L,Q(2:end-1,end),dx:dx:L,M(2:end-1,end)./(Q(2:end-1,end)-B(2:end-1)),dx:dx:L,B(2:end-1))
legend('height','M','B')
cons = sum(Q(:,:))/(xSteps+1);

end

