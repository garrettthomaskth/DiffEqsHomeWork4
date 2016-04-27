function [ phiLambda ] = phiDelta( lambda, delta )
%Phi Delta function for harten's hentropy fix
if abs(lambda) >= delta
    phiLambda = abs(lambda);
else
    phiLambda = (lambda^2 + delta^2)/(2*delta);
end

end

