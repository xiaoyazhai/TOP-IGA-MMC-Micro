function H=Heaviside(phi,alpha,epsilon)
H = zeros(size(phi));
H(phi>epsilon)=1;
H(phi<-epsilon)=alpha;
H(phi>=-epsilon & phi <= epsilon)=3*(1-alpha)/4.*(phi(phi>=-epsilon & phi <= epsilon)./epsilon - phi(phi>=-epsilon & phi <= epsilon).^3/(3*(epsilon)^3))+(1+alpha)/2;
end