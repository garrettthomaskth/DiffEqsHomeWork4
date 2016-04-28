function [ Ft ] = Ftilda( Lambda, W, j)

s = 0.5*(Lambda(j,1)+Lambda(j,2));
if Lambda(j,1) > 0 
    theta1 = W(j-1,1:2)*W(j,1:2)'/(W(j,1:2)*W(j,1:2)');
else
    theta1 = W(j+1,1:2)*W(j,1:2)'/(W(j,1:2)*W(j,1:2)');
end
W1 = W(j,1:2)*minmod(1,theta1);
fir = 0.5*abs(s)*(1-abs(s))*W1;

if Lambda(j,1) > 0 
    theta2 = W(j-1,3:4)*W(j,3:4)'/(W(j,3:4)*W(j,3:4)');
else
    theta2 = W(j+1,3:4)*W(j,3:4)'/(W(j,3:4)*W(j,3:4)');
end
W2 = W(j,1:2)*minmod(1,theta2);
sec = 0.5*abs(s)*(1-abs(s))*W2;

Ft = fir+sec;
end

