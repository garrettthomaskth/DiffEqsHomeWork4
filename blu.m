[Q,x,t,cons] = roeFirstFlat(200,600);
t = 6*[0 10 20 60 70 80];
for i=1:6
subplot(2,3,i);
plot(x,Q(:,t(i)+1));str = sprintf('T = %f',2.5/600*t(i));title(str);%axis([0 10 1 1.2]);
xlabel('x');ylabel('Height   h');
end