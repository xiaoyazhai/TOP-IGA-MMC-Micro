function [Ke,IndexCP,RBasisGP,numELE] = IGA_KE(D,px,py,eledof,KnotU1,KnotV1,U_KV,V_KV,CPxy,numCPX,numCPY,weights)
%--------------------------------------------------------------------------
% FUNCTION: Implement isogeometric analysis
% INPUT:
% D              : The constitutive matrix
% px,py          : Degree of basis functions
% eledof         : Degree of freedom per element
% KnotU1,KnotV1  : Non-repetitive knots in u/v-direction
% U_KU,V_KV      : Knot vectors in u/v-direction
% CPxy           : Control points
% numCPX,numCPY  : Number of control points in x/y-direction
% weights        : Weights of NURBS basis functions
% OUTPUT:
% Ke             : Element stiffness matrix
% IndexCP        : Indexes of control points
% RBasisGP       : NURBS basis functions at Gaussian points
% numELE         : Number of NURBS elements
%--------------------------------------------------------------------------
nGsx = 3; nGsy = 3; [GsIU,GsWU] = gspts(KnotU1,nGsx); [GsIV,GsWV] = gspts(KnotV1,nGsy); 
GsIW = reshape(kron(GsWU,GsWV),nGsx,nGsy,[]);
MidKnotU = (KnotU1(1:end-1) + KnotU1(2:end))/2; 
MidKnotV = (KnotV1(1:end-1) + KnotV1(2:end))/2;
IndexMKU = FindSpan(numCPX-1,px,MidKnotU,U_KV);
IndexMKV = FindSpan(numCPY-1,py,MidKnotV,V_KV);
IIndexMKU = repmat(IndexMKU,1,length(IndexMKV)); 
IIndexMKV = kron(IndexMKV,ones(1,length(IndexMKU))); 
numIndexCP = 0;
for icpy = py:-1:0
    for icpx = px:-1:0
        numIndexCP = numIndexCP + 1;
        IndexCP(numIndexCP,:) = (IIndexMKU - icpx)  + (IIndexMKV - icpy - 1) * numCPX;
    end
end
[numCPele,numELE]= size(IndexCP); CPWts = weights(IndexCP);
Ke = zeros(eledof,eledof,numELE); 
RBasisGP = zeros(nGsx*nGsy,numCPele,numELE);
for igsy = 1:nGsy
    derGsv = Der1BasisFuns(IndexMKV,GsIV(igsy,:),py,V_KV);
    NGsShapev = derGsv(1:py+1,:); DerGsShapev = derGsv(py+2:2*py+2,:);
    for igsx = 1:nGsx
        derGsu = Der1BasisFuns(IndexMKU,GsIU(igsx,:),px,U_KV);
        NGsShapeu = derGsu(1:px+1,:); DerGsShapeu = derGsu(px+2:2*px+2,:);
        WtGsP = kron(NGsShapev,NGsShapeu) .* CPWts;  sumWtGsP = sum(WtGsP); 
        WtGsPderu = kron(NGsShapev,DerGsShapeu) .* CPWts;  sumWtGsPderu = sum(WtGsPderu);
        WtGsPderv = kron(DerGsShapev,NGsShapeu) .* CPWts;  sumWtGsPderv = sum(WtGsPderv);
        % NURBS basis functions and their first partial derivatives
        RBasisGsP = WtGsP ./ repmat(sumWtGsP,(px+1)*(py+1),1);
        RBasisGsPderu = (WtGsPderu .* repmat(sumWtGsP,(px+1)*(py+1),1) - repmat(sumWtGsPderu,(px+1)*(py+1),1).* WtGsP)./ repmat(sumWtGsP.*sumWtGsP,(px+1)*(py+1),1);
        RBasisGsPderv = (WtGsPderv .* repmat(sumWtGsP,(px+1)*(py+1),1) - repmat(sumWtGsPderv,(px+1)*(py+1),1).* WtGsP)./ repmat(sumWtGsP.*sumWtGsP,(px+1)*(py+1),1);
        RBasisGsPDeruv(1,:,:) = RBasisGsPderu;  RBasisGsPDeruv(2,:,:) = RBasisGsPderv;
        RBasisGP((igsy - 1) * nGsx + igsx,:,:)=RBasisGsP; 
        % Control points of all NURBS elements 
        ELECP = zeros(numCPele,2,numELE);
        ELECP(:,1,:) = CPxy(IndexCP); ELECP(:,2,:) = CPxy(length(CPxy) + IndexCP);
        % Jacobian
        DerNxy = zeros(2,numCPele,numELE);
        eleGsPdetJa = zeros(numELE,1);
        for iele = 1:numELE
            eleGsPJcb = RBasisGsPDeruv(:,:,iele) * ELECP(:,:,iele); 
            eleGsPdetJa(iele) = det(eleGsPJcb); 
            DerNxy(:,:,iele) = eleGsPJcb\RBasisGsPDeruv(:,:,iele); 
        end
        % Strain-displacement matrix B 
        B = zeros(3,eledof,numELE);
        B(1,1:2:end,:) = DerNxy(1,:,:);
        B(2,2:2:end,:) = DerNxy(2,:,:);
        B(3,1:2:end,:) = DerNxy(2,:,:);
        B(3,2:2:end,:) = DerNxy(1,:,:);
        % Element stiffness matrix
        for ike = 1:numELE
            Ke(:,:,ike) = Ke(:,:,ike) + B(:,:,ike)'*D*B(:,:,ike)*eleGsPdetJa(ike)*GsIW(igsx,igsy,ike);
        end
    end
end