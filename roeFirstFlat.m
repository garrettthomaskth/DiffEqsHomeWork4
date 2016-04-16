function [ ] = roeFirstFlat(xSteps, tSteps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Chosen Variable
L=10; %L long
T=10;
%Given Variables
H=1;
w=0.1*L;
alpha=H/5;

dx=L/n;
dt=T/n;
%Initialize Matrix
Q = zeros(xSteps+2, 2*(tSteps+1));
%Initial Conditions
Q(2:(end-1),1) = H+epsilon*exp(-((dx/2:dx:(L-dx/2)) - L/2).^2/w^2);
%Choose initial conditions for momentum


f = @(u) [ u(2) , u(2)^2./u(1) + 0.5*g*u(1).^2]
F = @(u1,u2) 0.5*(f(u1)+f(u2)) - 0.5(
for i = 1:tSteps+1
    % Ghost point values
    Q(1,2*i-1) = Q(2,2*i-1);
    Q(1,2*i) = -Q(2,2*i);
    Q(end,2*i-1) = Q(end-1,2*i-1);
    Q(end,2*i) = -Q(end-1,2*i);
    
    for j = 1:xSteps+1
        F(j,:) = FLxF( Q(j+1,(2*i-1):(2*i)), Q(j,(2*i-1):(2*i)) );
    end
    for j = 2:xSteps+1
        Q(j,2*i+1:2*i+2) = Q(j,2*i-1:2*i) - dt/dx * (F(j,:)-F(j-1,:));
    end
end
x = linspace(0,L,xSteps);
t = linspace(0,Tend,tSteps+1);
Q = Q(2:end-1,2*(1:tSteps+1)-1);
cons = sum(Q(:,:))/(xSteps+1);
end

