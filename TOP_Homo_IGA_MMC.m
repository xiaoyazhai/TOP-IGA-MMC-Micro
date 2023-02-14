function TOP_Homo_IGA_MMC(lx,ly,nelx,nely,x_int,y_int,ini_val,volfrac,penal)
%--------------------------------------------------------------------------
% FUNCTION: Material design using the IGA-SIMP method
% INPUT:
% lx,ly     : Width and height of the base cell
% nelx,nely : Number of NURBS elements in x/y-direction
% volfrac   : The given upper limit of solid materials
% penal     : The penal factor
%--------------------------------------------------------------------------
close all
folder = [num2str(nelx),'_',num2str(volfrac),'_',num2str(penal),'_',num2str(ini_val(2)),'_',num2str(ini_val(5))];
mkdir(folder)

%% MATERIAL PROPERTIESs
CellV = lx*ly; 
% E0 = 71.7; nu = 0.33;
E0 = 1; nu = 0.3;
D = E0/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2.0];
%% PREPARE ISOGEOMETRIC ANALYSIS
% Initialize using NURBS patch
px = 2; py = 2; eledof = 2*(px+1)*(py+1);
numCPX = nelx + px; numCPY = nely + py;
U_Open = [zeros(1,px),linspace(0,1,nelx+1),ones(1,px)]; 
V_Open = [zeros(1,py),linspace(0,1,nely+1),ones(1,py)];
KnotU1 = unique(U_Open); KnotV1 = unique(V_Open);
IntervalU = KnotU1(2:end) - KnotU1(1:end-1);
IntervalV = KnotV1(2:end) - KnotV1(1:end-1);
CP_x0 = zeros(numCPX,1); CP_y0 = zeros(numCPY,1);
for gaxi = 1:(numCPX)
    CP_x0(gaxi) = sum(U_Open(gaxi+1:gaxi+px)/px)*lx;
end
for gayi = 1:(numCPY)
    CP_y0(gayi) = sum(V_Open(gayi+1:gayi+py)/py)*ly;
end
CPx = repmat(CP_x0,numCPY,1); CPy = kron(CP_y0,ones(numCPX,1));
CPgrid.x = reshape(CPx,numCPX,numCPY);
CPgrid.y = reshape(CPy,numCPX,numCPY);
CPxy = [CPx,CPy];
weights = ones(numCPX*numCPY,1);
% Component Geometry Initialization

% x0 = x_int/2:x_int:lx;
% y0 = y_int/2:y_int:ly;
% xn = length(x0);
% yn = length(y0);
% x0 = repmat(kron(x0,ones(1,2)),1,yn);
% y0 = kron(y0,ones(1,2*xn));

tempHx0 = x_int/2:x_int:lx;
tempHy0 = 0:ini_val(1):lx;
% 2*2
% x0 = [tempHx0,tempHy0,tempHx0,tempHy0,tempHx0]; %,tempHy0,tempHx0,tempHy0,tempHx0
% y0 = [repmat(tempHy0(1),1,2),repmat(tempHx0(1),1,3),repmat(tempHy0(2),1,2),...
%     repmat(tempHx0(2),1,3),repmat(tempHy0(3),1,2)];

% 3*3
x0 = [tempHx0,tempHy0,tempHx0,tempHy0,tempHx0,tempHy0,tempHx0]; 
y0 = [repmat(tempHy0(1),1,3),repmat(tempHx0(1),1,4),repmat(tempHy0(2),1,3),...
    repmat(tempHx0(2),1,4),repmat(tempHy0(3),1,3),repmat(tempHx0(3),1,4),repmat(tempHy0(4),1,3)];

% 4*4 
% x0 = [tempHx0,tempHy0,tempHx0,tempHy0,tempHx0,tempHy0,tempHx0,tempHy0,tempHx0];
% y0 = [repmat(tempHy0(1),1,4),repmat(tempHx0(1),1,5),repmat(tempHy0(2),1,4),...
%     repmat(tempHx0(2),1,5),repmat(tempHy0(3),1,4),repmat(tempHx0(3),1,5),...
%     repmat(tempHy0(4),1,4),repmat(tempHx0(4),1,5),repmat(tempHy0(5),1,4)];

NC = length(x0);
L = repmat(ini_val(1),1,NC);
t1 = repmat(ini_val(2),1,NC);
t2 = repmat(ini_val(3),1,NC);
t3 = repmat(ini_val(4),1,NC);

% st = repmat([ini_val(5) -ini_val(5)],1,NC/2);
% 2*2
% st = [repmat(ini_val(5),1,2),repmat(sqrt(1-ini_val(5)),1,3),repmat(ini_val(5),1,2),repmat(sqrt(1-ini_val(5)),1,3),repmat(ini_val(5),1,2)];

% 3*3
st = [repmat(ini_val(5),1,3),repmat(sqrt(1-ini_val(5)),1,4),repmat(ini_val(5),1,3),...
    repmat(sqrt(1-ini_val(5)),1,4),repmat(ini_val(5),1,3),repmat(sqrt(1-ini_val(5)),1,4),repmat(ini_val(5),1,3)];

% 4*4
% st = [repmat(ini_val(5),1,4),repmat(sqrt(1-ini_val(5)),1,5),repmat(ini_val(5),1,4),repmat(sqrt(1-ini_val(5)),1,5),repmat(ini_val(5),1,4),...
%     repmat(sqrt(1-ini_val(5)),1,5),repmat(ini_val(5),1,4),repmat(sqrt(1-ini_val(5)),1,5),repmat(ini_val(5),1,4)];

variable = [x0;y0;L;t1;t2;t3;st];
% Parameter of MMA
xval = variable(:);
xold1 = xval;
xold2 = xval;
% Limits of variables
xmin = [0; 0; 0.001; 0.001; 0.001; 0.03; -1.0];
xmin = repmat(xmin,NC,1);
xmax = [lx; ly; 2.0; 0.2; 0.2; 0.2; 1.0];
xmax = repmat(xmax,NC,1);
low = xmin;
upp = xmax;
m = 1;
Var_num = 7;
nn = Var_num*NC;
c = 1000*ones(m,1);
d = zeros(m,1);
a0 = 1;
a = zeros(m,1);
% Stiffness matrix
[KE,IndexCP,RBasisGP,numELE] = IGA_KE(D,px,py,eledof,KnotU1,KnotV1,U_Open,V_Open,CPxy,numCPX,numCPY,weights);
IndexCPdof(1:2:eledof,:) = 2 * IndexCP - 1; IndexCPdof(2:2:eledof,:) = 2 * IndexCP; 
iK = kron(IndexCPdof,ones(1,eledof));
jK = kron(IndexCPdof,ones(eledof,1));
%% PERIODIC BOUNDARY CONDITIONS
epsilon = eye(3); ufixedCP = zeros(8,3);
nICP = numCPX*numCPY; ndof = 2*nICP;
alldofs = 1:ndof; chi = zeros(ndof,3);
nCPArray = reshape(1:numCPX*numCPY,numCPX,numCPY);
n1 = [nCPArray([1,end],1);nCPArray([1,end],end)];    %Four corners
d1 = reshape([(2*n1-1),2*n1]',1,8);
n3 = [nCPArray(2:end-1,1)',nCPArray(1,2:end-1)];     %Bottom and left edges
d3 = reshape([2*n3-1;2*n3],1,2*(numCPX+numCPY-4));
n4 = [nCPArray(2:end-1,end)',nCPArray(end,2:end-1)]; %Top and right edges
d4 = reshape([2*n4-1;2*n4],1,2*(numCPX+numCPY-4));
d2 = setdiff(alldofs,[d1,d3,d4]); %Inner control points
for ie = 1:3
    Eepsilon = [epsilon(1,ie), epsilon(3,ie)/2; epsilon(3,ie)/2, epsilon(2,ie)];
    ufixedCP(3:4,ie) = Eepsilon * [lx; 0];
    ufixedCP(5:6,ie) = Eepsilon * [0; ly];
    ufixedCP(7:8,ie) = ufixedCP(3:4,ie) + ufixedCP(5:6,ie);
end
wfixedCP = [repmat(ufixedCP(5:6,:),numCPX-2,1); repmat(ufixedCP(3:4,:),numCPY-2,1)]; 
%% Iteration Initialization
p = 6;
alpha = 1e-3;
epsilon = 0.1;
Phi = zeros(numCPX,numCPY,NC);
Loop = 0;
Change = 1.0;
MaxIter = 500;
MinIter = 50;
IterComp = zeros(1,MaxIter);
IterVolume = zeros(1,MaxIter);
%% Start Iteration
tic
while (Change > 0.001 && Loop < MaxIter) || Loop < MinIter
    Loop = Loop + 1;
    for ni = 1:NC
        Phi(:,:,ni) = tPhi(variable(:,ni),CPgrid.x,CPgrid.y,p);
    end
    Phi_max = max(Phi,[],3);
    % Plot Components
    h1 = figure(1);
    clf(h1)
   % hold on
    contourf(CPgrid.x,CPgrid.y,Phi_max,[0,0]);
    %contourf(CPgrid.x+1,CPgrid.y,Phi_max,[0,0]);
    %contourf(CPgrid.x,CPgrid.y+1,Phi_max,[0,0]);
    %contourf(CPgrid.x+1,CPgrid.y+1,Phi_max,[0,0]);
   % hold off
    axis equal; axis([0 lx 0 ly]); 
    xticks(0:0.5:lx);
    yticks(0:0.5:ly);
    if Loop == 1
        print(h1,[folder,'\44MMC-1-',num2str(lx),'-',num2str(ly),'-',num2str(nelx),'-',num2str(nely),'-',num2str(volfrac),'-',...
            num2str(ini_val(2)),'-',num2str(ini_val(3)),'-',num2str(ini_val(4)),'(',num2str(1),')','.png'],'-dpng');
    end
    % Calculating the Finite Difference Quotient of Heaviside Function H
    H = Heaviside(Phi_max,alpha,epsilon);
    diffH = cell(NC,1);
    deltaHx = max(2*min(min(IntervalU),min(IntervalV)),0.005);%做有限差分时，那个改变的小量
    for j = 1:NC
        for ii = 1:Var_num
            variable1 = variable;
            tempPhi1 = Phi;
            variable1(ii,j) = variable(ii,j) + deltaHx;
            tempPhiD1 = tPhi(variable1(:,j),CPgrid.x,CPgrid.y,p);
            tempPhi1(:,:,j) = tempPhiD1;
            tempPhi_max1 = max(tempPhi1,[],3);
            
            variable2 = variable;
            tempPhi2 = Phi;
            variable2(ii,j) = variable(ii,j) - deltaHx;
            tempPhiD2 = tPhi(variable2(:,j),CPgrid.x,CPgrid.y,p);
            tempPhi2(:,:,j) = tempPhiD2;
            tempPhi_max2 = max(tempPhi2,[],3);
            
            HD1 = Heaviside(tempPhi_max1,alpha,epsilon);
            HD2 = Heaviside(tempPhi_max2,alpha,epsilon);
            tempdiffH= (HD1-HD2)/(2*deltaHx);
            diffH{j}(:,ii) = tempdiffH(:);
        end
    end
    
    TempdenkModulus=sum(squeeze(sum(RBasisGP.*reshape(kron(H(IndexCP),ones(1,(px+1)*(py+1))),(px+1)*(py+1),(px+1)*(py+1),[]))))/size(RBasisGP,1); %这里TempdenkModulus表示：用NURBS基函数对弹性模量进行插值，得到每个高斯点处的弹性模量，并进行均分，得到单元对应的弹性模量插值系数
    denk = TempdenkModulus.^penal ; %带惩罚的单元刚度阵插值系数
    krondenk = kron(denk,ones(eledof,eledof));
    threeDdenk = reshape(krondenk,eledof,eledof,[]);
    sK = KE.*threeDdenk;
    K = sparse(iK(:),jK(:),sK(:)); K = (K+K')/2;
    Kr = [K(d2,d2),K(d2,d3)+K(d2,d4); K(d3,d2)+K(d4,d2), K(d3,d3)+K(d4,d3)+K(d3,d4)+K(d4,d4)];
    Fr = -[K(d2,d1); K(d3,d1)+K(d4,d1)] * ufixedCP - [K(d2,d4); K(d3,d4)+K(d4,d4)] * wfixedCP;
    chi(d1,:) = ufixedCP;
    chi([d2,d3],:) = Kr\Fr;
    chi(d4,:) = chi(d3,:) + wfixedCP;
    ChiE = zeros(3,3,numELE);
    for iele = 1:numELE
        tempchi = chi(IndexCPdof(:,iele),:);
        ChiE(:,:,iele) = tempchi' * sK(:,:,iele) * tempchi;
    end
    CHI = sum(ChiE,3)/CellV;
    %% Objective function and sensitivity analysis
    Comp =  CHI(1,2) - (0.8^Loop) * (CHI(1,1)+CHI(2,2)); % -CHI(3,3);%
    
    % Volume 
    den = sum(H(IndexCP)) / ((px+1)*(py+1)); %体积插值系数
    repmatIntervalU = repmat(IntervalU,1,length(KnotV1)-1)*lx;
    kronIntervalV = kron(IntervalV,ones(1,length(KnotU1)-1))*ly;
    VolumeA1 = sum(den.*repmatIntervalU.*kronIntervalV);
    % Sensitivity analysis
    df0dx = zeros(NC * Var_num,1);
    dfdx = zeros(NC * Var_num,1);
    % Energy_element %这里的单元能量需要重新求
    % Energy_element = sum(squeeze(sum(reshape(kron(U_Open(IndexCPdof),ones(1,eledof)),eledof,eledof,[]).*KE)).*U(IndexCPdof));%计算每个单元对应的应变能 %在这里把厚度放在的KE内 %未插值处的应变能
    for j = 1:NC
        for i = 1:Var_num
            jCompiVarH = diffH{j}(:,i); %第j个组件，第i个设计变量，对应于每个高斯积分点的diffH值
            % tempParDer1 = sum(RBasisGP.*reshape(kron(jCompiVarH(IndexCP),ones(1,(px+1)*(py+1))),(px+1)*(py+1),(px+1)*(py+1),[]));
            ParDer = sum(squeeze(sum(RBasisGP.*reshape(kron(jCompiVarH(IndexCP),ones(1,(px+1)*(py+1))),(px+1)*(py+1),(px+1)*(py+1),[]))))/size(RBasisGP,1);% .*WeightELEGP(:));
            derKcoe1 = penal*(TempdenkModulus.^(penal-1)).*ParDer;
            tempderKcoe = kron(derKcoe1,ones(3,3)); derKcoe3 = reshape(tempderKcoe,3,3,[]);
            tempdf0dx = derKcoe3 .* ChiE; 
            dIcCHI = sum(tempdf0dx,3)/CellV;
            df0dx((j-1)*Var_num + i) = dIcCHI(1,2) - (0.8^Loop) * (dIcCHI(1,1)+dIcCHI(2,2)); % -dIcCHI(3,3);%
            dfdx((j-1)*Var_num + i) = sum(sum(jCompiVarH(IndexCP))/((px+1)*(py+1)*nelx*nely));
        end
    end
    
    % MMA Optimization
    f0val = Comp;
    df0dx= df0dx/max(abs(df0dx));
    fval = VolumeA1/(lx*ly)-volfrac;
    dfdx= dfdx/max(abs(dfdx));
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,ss,low,upp] = ...
        mmasub(m,nn,Loop,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
    xold2 = xold1;
    xold1 = xval;
    Change=max(abs(xval-xmma));
    if Loop < MinIter
        xval = xmma;
    elseif Loop < MaxIter/2
        xval = 0.4*xmma + 0.6*xold1;
    else
        xval = 0.3*xmma + 0.6*xold1 + 0.1*xold2;
    end
    variable = reshape(round(xval*1e4)/1e4,Var_num,NC);
    disp([' It.: ' sprintf('%4i\t',Loop) ' Obj.: ' sprintf('%6.4f\t',f0val) ' Nu.: ' sprintf('%6.4f\t',CHI(1,2)/CHI(1,1)) ' Vol.: ' ...
        sprintf('%6.4f\t',fval) 'Ch.:' sprintf('%6.4f\t',Change)]);
    xlswrite([folder,'\',folder,'.xlsx'],[Loop,f0val,CHI(1,2)/CHI(1,1),fval,Change,toc],1,['A',num2str(Loop)]);
    if(mod(Loop,3)==1)
        saveas(1,[folder,'\',folder,'_',num2str(Loop),'.png']);
    end
end
Loop = Loop + 1;
for ni = 1:NC
    Phi(:,:,ni) = tPhi(variable(:,ni),CPgrid.x,CPgrid.y,p);
end
Phi_max = max(Phi,[],3);
h1 = figure(1);
contourf(CPgrid.x,CPgrid.y,Phi_max,[0,0]);
axis equal; axis([0 lx 0 ly]);
xticks(0:0.5:lx);
yticks(0:0.5:ly);
print(h1,[folder,'\',num2str(fval),'_',num2str(CHI(1,2)/CHI(1,1)),'.png'],'-dpng');

H = Heaviside(Phi_max,alpha,epsilon);
TempdenkModulus=sum(squeeze(sum(RBasisGP.*reshape(kron(H(IndexCP),ones(1,(px+1)*(py+1))),(px+1)*(py+1),(px+1)*(py+1),[]))))/size(RBasisGP,1); %这里TempdenkModulus表示：用NURBS基函数对弹性模量进行插值，得到每个高斯点处的弹性模量，并进行均分，得到单元对应的弹性模量插值系数
denk = TempdenkModulus.^penal ; %带惩罚的单元刚度阵插值系数
krondenk = kron(denk,ones(eledof,eledof));
threeDdenk = reshape(krondenk,eledof,eledof,[]);
sK = KE.*threeDdenk;
K = sparse(iK(:),jK(:),sK(:)); K = (K+K')/2;
Kr = [K(d2,d2),K(d2,d3)+K(d2,d4); K(d3,d2)+K(d4,d2), K(d3,d3)+K(d4,d3)+K(d3,d4)+K(d4,d4)];
Fr = -[K(d2,d1); K(d3,d1)+K(d4,d1)] * ufixedCP - [K(d2,d4); K(d3,d4)+K(d4,d4)] * wfixedCP;
chi(d1,:) = ufixedCP;
chi([d2,d3],:) = Kr\Fr;
chi(d4,:) = chi(d3,:) + wfixedCP;
ChiE = zeros(3,3,numELE);
for iele = 1:numELE
    tempchi = chi(IndexCPdof(:,iele),:);
    ChiE(:,:,iele) = tempchi' * sK(:,:,iele) * tempchi;
end
CHI = sum(ChiE,3)/CellV;
%% Objective function and sensitivity analysis
Comp =  CHI(1,2) - (0.8^Loop) * (CHI(1,1)+CHI(2,2)); % -CHI(3,3);%

% Volume
den = sum(H(IndexCP)) / ((px+1)*(py+1)); %体积插值系数
repmatIntervalU = repmat(IntervalU,1,length(KnotV1)-1)*lx;
kronIntervalV = kron(IntervalV,ones(1,length(KnotU1)-1))*ly;
VolumeA1 = sum(den.*repmatIntervalU.*kronIntervalV);
f0val = Comp;
fval = VolumeA1/(lx*ly)-volfrac;
disp([' It.: ' sprintf('%4i\t',Loop) ' Obj.: ' sprintf('%6.4f\t',f0val) ' Nu.: ' sprintf('%6.4f\t',CHI(1,2)/CHI(1,1)) ' Vol.: ' ...
        sprintf('%6.4f\t',fval) 'Ch.:' sprintf('%6.4f\t',Change)]);

h2= figure(2);
hold on
contourf(CPgrid.x,CPgrid.y,Phi_max,[0,0]);
contourf(CPgrid.x+1,CPgrid.y,Phi_max,[0,0]);
contourf(CPgrid.x,CPgrid.y+1,Phi_max,[0,0]);
contourf(CPgrid.x+1,CPgrid.y+1,Phi_max,[0,0]);
hold off
axis equal; axis([0 2*lx 0 2*ly]);
xticks(0:0.5:2*lx);
yticks(0:0.5:2*ly);
print(h2,[folder,'/44MMC-4-',num2str(lx),'-',num2str(ly),'-',num2str(nelx),'-',num2str(nely),'-',num2str(volfrac),'-',...
            num2str(ini_val(2)),'-',num2str(ini_val(3)),'-',num2str(ini_val(4)),'(',num2str(1),')','.png'],'-dpng');
        
h3= figure(3);
hold on
contourf(CPgrid.x,CPgrid.y,Phi_max,[0,0]);
contourf(CPgrid.x+1,CPgrid.y,Phi_max,[0,0]);
contourf(CPgrid.x+2,CPgrid.y,Phi_max,[0,0]);
contourf(CPgrid.x,CPgrid.y+1,Phi_max,[0,0]);
contourf(CPgrid.x+1,CPgrid.y+1,Phi_max,[0,0]);
contourf(CPgrid.x+2,CPgrid.y+1,Phi_max,[0,0]);
contourf(CPgrid.x,CPgrid.y+2,Phi_max,[0,0]);
contourf(CPgrid.x+1,CPgrid.y+2,Phi_max,[0,0]);
contourf(CPgrid.x+2,CPgrid.y+2,Phi_max,[0,0]);
hold off
axis equal; axis([0 3*lx 0 3*ly]);
xticks(0:0.5:3*lx);
yticks(0:0.5:3*ly);
print(h3,[folder,'/44MMC-9-',num2str(lx),'-',num2str(ly),'-',num2str(nelx),'-',num2str(nely),'-',num2str(volfrac),'-',...
            num2str(ini_val(2)),'-',num2str(ini_val(3)),'-',num2str(ini_val(4)),'(',num2str(1),')','.png'],'-dpng');


