function [Q,U,x,t,cons] = roeNotFlat(xSteps, tSteps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Chosen Variable
L=10; %L long
T=15;


%Given Variables
H=1;
w=0.1*L;
alpha=H/5;
g=9.8;

r=L/6;
B0 = H/10;

dx=L/xSteps;
dt=T/tSteps;
%Initialize Matrix
Q = zeros(xSteps+2, 2*(tSteps+1));
%Initial Conditions
Q(1:(end),1) = 0.8;
Q(1:(end),2) = 0.88;
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
%Q(:,2) = -(Q(:,1)-H).*sqrt(g*(Q(:,1)));


f = @(u) [ u(2) , u(2)^2./u(1) + 0.5*g*u(1).^2];
%Ffun = @(u1,u2) 0.5*(f(u1)+f(u2)) - 0.5*(abs(lambda1)*W1 + abs(lambda2)*W2);

F = zeros(xSteps+1,2);
S = zeros(xSteps+2,2);
for i = 1:tSteps+1
    % Ghost point values
    Q(1,2*i-1) = 0.8;
    Q(1,2*i) = 0.88;
    Q(end,2*i-1:2*i) = 2*Q(end-1,2*i-1:2*i) - Q(end-2,2*i-1:2*i);
    
    for j = 1:xSteps+1
        cHat = sqrt(0.5*g*(Q(j,2*i-1)+Q(j+1,2*i-1)));
        uHat = (sqrt(Q(j,2*i-1))*Q(j,2*i)/Q(j,2*i-1)+sqrt(Q(j+1,2*i-1))*Q(j+1,2*i)/Q(j+1,2*i-1))/(sqrt(Q(j,2*i-1))+sqrt(Q(j+1,2*i-1)));
        
        delta = Q(j+1,2*i-1:2*i)-Q(j,2*i-1:2*i);
        
        W1 = ((uHat+cHat)*delta(1)-delta(2))/(2*cHat) * [1, uHat - cHat];
        W2 = (-(uHat-cHat)*delta(1)+delta(2))/(2*cHat) * [1, uHat + cHat];
        
        lambda1 = uHat - cHat;
        lambda2 = uHat + cHat;
        
        F(j,:) = 0.5*(f(Q(j+1,(2*i-1):(2*i)))+f(Q(j,(2*i-1):(2*i)))) - 0.5*(abs(lambda1)*W1 + abs(lambda2)*W2);
        
        x = (j-1.5)*dx;
        
        if abs(x - L/2) < r            
            S(j,2) = g*(Q(j,2*i-1))*pi*B0*sin(pi*(x-L/2)/r)/(2*r);
        end
    end
    for j = 2:xSteps+1
        Q(j,2*i+1:2*i+2) = Q(j,2*i-1:2*i) - dt/dx * (F(j,:)-F(j-1,:))+ dt*S(j,:);
    end
end
for j = 2:xSteps+1
    x = (j-1.5)*dx;
        if abs(x - L/2) < r            
             Q(j,2*(1:tSteps+1)-1) = Q(j,2*(1:tSteps+1)-1)+B0*cos(pi*(x-L/2)/(2*r))^2;
        end
end
x = dx/2:dx:L-dx/2;
t = linspace(0,T,tSteps+1);
U = Q(2:end-1,2*(1:tSteps+1));
Q = Q(2:end-1,2*(1:tSteps+1)-1);
mesh(t,x,Q)
rotate3d on
cons = sum(Q(:,:))/(xSteps+1);
end

