


[Q1,x1,t,cons] = roeFirstFlat(80,80);
[Q2,x2,t,cons] = roeFirstFlat(160,160);
[Q3,x3,t,cons] = roeFirstFlat(320,320);

t = [0 8 16 24 36 48 60 72 80];
for i=1:9
subplot(3,3,i);
plot(x1,Q1(:,t(i)+1),x2,Q2(:,2*t(i)+1),x3,Q3(:,4*t(i)+1));str = sprintf('T = %f',2.5/80*t(i));title(str);legend('N=80','N=160','N=320');
xlabel('x');ylabel('Height   h');
end