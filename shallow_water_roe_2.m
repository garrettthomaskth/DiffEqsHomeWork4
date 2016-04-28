function [  ] = shallow_water_roe_2(  )
% Shallow water model solved with Roe scheme

C=0.2; % Stability parameter
N=100; % Number of space cells
s=600; % Number of time steps
h=zeros(N,s);
m=zeros(N,s);
g=9.61;


% Definie parameters
H=1; L=10; w=0.1*L; a=H/5; %a = H/5;
dx=L/N;
dt= C*dx;

% Bottom (=0 for horizontal)
Bx= zeros(N,1);
BB=zeros(N,1);
B0=H/10;
r=L/6;
jmin = ceil((L/2-r)/dx);
jmax = floor((L/2+r)/dx);
for j= jmin:jmax
    Bx(j)= -B0*pi/r*sin(pi/(2*r)*((j-1/2)*dx-L/2))*cos(pi/(2*r)*((j-1/2)*dx-L/2));
    BB(j) = B0*(cos(pi*((j-1/2)*dx-L/2)/(2*r)))^2;
end

BB
for j=1:N %Initial condition
   h(j,1)= H+a*exp(-(dx*(j-1/2)-L/2)^2/w^2);
   %m(j,1) = h(j,1)^(3/2)*sqrt(g);
   m(j,1) = (h(j,1)-H+BB(j))*sqrt(g*(h(j,1)))+0.5;
end



% Solve every time step
for n=1:s-1
    
    % Add BC
    %he=[h(1,n); h(:,n); h(N,n)];
    %me=[-m(1,n); m(:,n); -m(N,n)];
    he=[1 ; h(:,n); 1];
    me=[m(1,n); m(:,n); m(N,n)];
    
    F=zeros(N+1,2);

    for j=1:N+1 
      f = [ me(j)+me(j+1), me(j)^2/he(j)+ 1/2*g*he(j)^2 + me(j+1)^2/he(j+1)+ 1/2*g*he(j+1)^2]; %f(Q_i)+f(Q_i+1)
       uhat= real((me(j)/sqrt(he(j))+ me(j+1)/sqrt(he(j+1)))/(sqrt(he(j))+sqrt(he(j+1))));
       chat= real(sqrt(g*(he(j+1)+he(j))/2));
       delta = [he(j+1)-he(j); me(j+1)-me(j)];
       alpha1 = ((uhat+chat)*delta(1)-delta(2))/(2*chat);
       alpha2 = (-(uhat-chat)*delta(1) + delta(2))/(2*chat);
       w1 = alpha1*[1, uhat-chat];
       w2 = alpha2*[1, uhat+chat];
    
      F(j,:) = 1/2*( f- abs(uhat-chat)*w1-abs(uhat+chat)*w2); %flux F_(j+1/2)
    end
    
    % forward Euler
    new= [he(2:N+1),me(2:N+1)] - dt/dx*(F(2:N+1,:)-F(1:N,:))+ dt*[zeros(N,1), -g*(he(2:N+1)).*Bx];
    
    % updating solution
    h(:,n+1) =  new(:,1);
    m(:,n+1) =  new(:,2); 
    
end


T=dt*s;
st=sprintf('Solution up to time T=%f', T)
dt

h=h+kron(BB, ones(1,s));

v= m./h;

% Plot solutions
t=linspace(0,T,s);
x=linspace(0,10,N);

figure(1)
mesh(t,x,h);
xlabel('time t')
ylabel('space x')
zlabel('height h')
title('Height of the water')
rotate3d on;

figure(2)
mesh(t,x,v)
xlabel('time t')
ylabel('space x')
zlabel('velocity v')
title('Velocity of the water')
rotate3d on;



for k=1:s 
    j=k;
    figure(3)
    plot(x,h(:,j));
    hold on;
    plot(x,BB);
    plot(x,v(:,j));
    xlabel('space x')
    ylabel('height')
    str=sprintf('Height of the water at time t=%f',(j-1)*dt);
    title(str)
    legend('Water height h+B', 'ground B', 'Velocity v','Location','West');
    hold off;
    %ylim([1,1.09])
    
    

    pause;
end

end