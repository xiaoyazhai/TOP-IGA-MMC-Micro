function ki = FindSpan(n,p,u,U)
%--------------------------------------------------------------------------
% FUNCTION: Find the knot span index for parameters u
% INPUT:
% n : Number of basis function -1
% u : Parameters to be evaluated
% U : Knot vector(row vector)
% OUTPUT:
% ki: Index of knot span
%--------------------------------------------------------------------------
if (max(u(:))>U(end) || min(u(:))<U(1))
  error('Some parameters are out of bounds');
end
ki = zeros(size(u));
u2 = find(u==U(n+2)); ki(u2) = n+1;
u1 = setdiff(1:length(u),u2); setdu = u(u1);
low = p*ones(size(setdu));
high = (n+1)*ones(size(setdu));
mid = floor((low+high)/2);
while max( setdu<U(mid+1) | setdu >= U(mid+2))==1
    high(setdu<U(mid+1)) = mid(setdu<U(mid+1));
    low(setdu>=U(mid+2)) = mid(setdu>=U(mid+2));
    mid = floor((low+high)/2);
end
ki(u1) = mid+1;