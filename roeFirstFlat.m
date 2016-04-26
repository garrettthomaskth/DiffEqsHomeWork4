function [Q,x,t,cons] = roeFirstFlat(xSteps, tSteps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Chosen Variable
L=10; %L long
T=10;
%Given Variables
H=1;
w=0.1*L;
alpha=H/5;
g=9.8;


dx=L/xSteps;
dt=T/tSteps;
%Initialize Matrix
Q = zeros(xSteps+2, 2*(tSteps+1));
%Initial Conditions
Q(2:(end-1),1) = H+alpha*exp(-((dx/2:dx:(L-dx/2)) - L/2).^2/w^2);
i=1;
Q(1,2*i-1) = Q(2,2*i-1);
Q(1,2*i) = -Q(2,2*i);
Q(end,2*i-1) = Q(end-1,2*i-1);
Q(end,2*i) = -Q(end-1,2*i);
%Choose initial conditions for momentum
Q(:,2) = -(Q(:,1)-H).*sqrt(g*(Q(:,1)));
%-(Q(:,1)-H).*sqrt(g*(Q(:,1)))


f = @(u) [ u(2) , u(2)^2./u(1) + 0.5*g*u(1).^2];
%Ffun = @(u1,u2) 0.5*(f(u1)+f(u2)) - 0.5*(abs(lambda1)*W1 + abs(lambda2)*W2);

F = zeros(xSteps+1,2);
for i = 1:tSteps+1
    % Ghost point values
    Q(1,2*i-1) = Q(2,2*i-1);
    Q(1,2*i) = -Q(2,2*i);
    Q(end,2*i-1) = Q(end-1,2*i-1);
    Q(end,2*i) = -Q(end-1,2*i);
    
    for j = 1:xSteps+1
        cHat = sqrt(0.5*g*(Q(j,2*i-1)+Q(j+1,2*i-1)));
        uHat = (sqrt(Q(j,2*i-1))*Q(j,2*i)/Q(j,2*i-1)+sqrt(Q(j+1,2*i-1))*Q(j+1,2*i)/Q(j+1,2*i-1))/(sqrt(Q(j,2*i-1))+sqrt(Q(j+1,2*i-1)));
        
        delta = Q(j+1,2*i-1:2*i)-Q(j,2*i-1:2*i);
        
        W1 = ((uHat+cHat)*delta(1)-delta(2))/(2*cHat) * [1, uHat - cHat];
        W2 = (-(uHat-cHat)*delta(1)+delta(2))/(2*cHat) * [1, uHat + cHat];
        
        lambda1 = uHat - cHat;
        lambda2 = uHat + cHat;
        
        F(j,:) = 0.5*(f(Q(j+1,(2*i-1):(2*i)))+f(Q(j,(2*i-1):(2*i)))) - 0.5*(abs(lambda1)*W1 + abs(lambda2)*W2);
        %Ffun( Q(j+1,(2*i-1):(2*i)), Q(j,(2*i-1):(2*i)) );
    end
    for j = 2:xSteps+1
        Q(j,2*i+1:2*i+2) = Q(j,2*i-1:2*i) - dt/dx * (F(j,:)-F(j-1,:));
    end
end
x = linspace(0,L,xSteps);
t = linspace(0,T,tSteps+1);
Q = Q(2:end-1,2*(1:tSteps+1)-1);
mesh(t,x,Q)
rotate3d on
cons = sum(Q(:,:))/(xSteps+1);
end

