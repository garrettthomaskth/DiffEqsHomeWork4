n=320;
[Q1,x1,t,cons] = advComp(n,0.25);
[Q2,x2,t,cons] = roeFirstFlat(n,n);
t = 4*[0 8 16 24 36 48 60 72 80];

for i=1:9
subplot(3,3,i);
plot(x1,Q1(:,t(i)+1),x2,Q2(:,t(i)+1));str = sprintf('T = %f',2.5/320*t(i));title(str);legend('Lax-Friedrichs','Roe');
xlabel('x');ylabel('Height   h');
end