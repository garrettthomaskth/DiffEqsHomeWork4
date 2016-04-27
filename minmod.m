function [ mm ] = minmod( a, b )
mm = a;
if abs(a) > abs(b) && a*b > 0
    mm = b;
elseif a*b <= 0
    mm = 0;
end

end

