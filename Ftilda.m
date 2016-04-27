function [ Ft ] = Ftilda( Lambda, W, j)

if Lambda(j,1) > 0 
    theta1 = W(j-1,1:2)*W(j,1:2)'/(W(j,1:2)*W(j,1:2)');
else
    theta1 = W(j+1,1:2)*W(j,1:2)'/(W(j,1:2)*W(j,1:2)');
end
W1 = W(j,1:2)*minmod(1,theta1);
fir = 0.05*abs(Lambda(j,1))*(1-abs(Lambda(j,1)))*W1;

if Lambda(j,1) > 0 
    theta2 = W(j-1,3:4)*W(j,3:4)'/(W(j,3:4)*W(j,3:4)');
else
    theta2 = W(j+1,3:4)*W(j,3:4)'/(W(j,3:4)*W(j,3:4)');
end
W2 = W(j,1:2)*minmod(1,theta2);
sec = 0.05*abs(Lambda(j,2))*(1-abs(Lambda(j,2)))*W2;

Ft = fir+sec;
end

