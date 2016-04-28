function [Q,x,t,cons] = adv(xSteps,ratio,Tend,alpha,epsilon)
%Problem 2.1 Lax-Friedrich method

%given variables
L = 10;
H = 1;
g = 9.61;
w = 0.4;
%calculate number of steps
dx = L/xSteps;
dt = ratio*dx;
tSteps = round(Tend/dt);
%flux function
f = @(u) [ u(2) , u(2)^2./u(1) + 0.5*g*u(1).^2];
%Lax-Friedrich flux function
FLxF = @(u2,u1) (0.5*(f(u2)+f(u1) - alpha*dx/dt*(u2-u1)));

%Initial Conditions
Q = zeros(xSteps+2, 2*(tSteps+1));
Q(2:(end-1),1) = H+epsilon*exp(-((dx/2:dx:(L-dx/2)) - L/2).^2/w^2);

F = zeros(xSteps+1,2);
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
x = linspace(0,L,xSteps+1);
t = linspace(0,Tend,tSteps+1);
Q = Q(1:end-1,2*(1:tSteps+1)-1);
cons = sum(Q(:,:))/(xSteps+1);
end
