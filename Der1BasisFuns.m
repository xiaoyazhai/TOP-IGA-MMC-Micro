function derfs = Der1BasisFuns(ki,u,p,U)
%--------------------------------------------------------------------------
% FUNCTION: Calculate B-spline basis functions and their first derivatives
% INPUT:
% ki   : corresponding knot indexes
% u    : parameters to be evaluated
% p    : degree of the B-spline basis functions
% U    : the knot vector
% OUTPUT:
% derfs: matrix (2*(p+1), length(u))
%--------------------------------------------------------------------------
numu = length(u);
derfs = zeros(2*(p+1),numu);
N = zeros(p+1,p+1,numu);
N(1,1,:) = 1;
left = zeros(p+1,numu);
right = zeros(p+1,numu);
for j = 1:p
   left(j+1,:) = u - U(ki-j+1); 
   right(j+1,:) =  U(ki+j) - u;
   saved = zeros(1,numu);
   for r = 0:j-1
       N(j+1,r+1,:) = right(r+2,:) + left(j-r+1,:);
       temp = reshape(N(r+1,j,:)./N(j+1,r+1,:),1,numu);
       N(r+1,j+1,:) = saved + right(r+2,:).*temp;
       saved = left(j-r+1,:).*temp;
   end
   N(j+1,j+1,:) = saved;
end
derfs(1:p+1,:) = reshape(N(:,p+1,:),p+1,numu);
for r = 0:p
   if (r>=1)
       derfs((p+1)+ (r+1),:) = reshape(N(r,p,:)./N(p+1,r,:),1,numu);
   end
   if (r<=p-1)
       derfs((p+1)+(r+1),:) = derfs(p+1+(r+1),:) - reshape(N(r+1,p,:)./N(p+1,r+1,:),1,numu);
   end
end
derfs(p+1+1:end,:) = derfs(p+1+1:end,:) * p;