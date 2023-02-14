function [Gs,Gw] = gspts(Vectors,numGaussP) 
%--------------------------------------------------------------------------
% FUNCTION: Calculate the Gaussian points and weights
% INPUT:
% Vectors   : Knot vectors (row vectors)
% numGaussP : Number of Gaussian points in an element
% OUTPUT:
% Gs        : Gaussian points
% Gw        : Gaussian weights
%--------------------------------------------------------------------------
h = (Vectors(2:end)-Vectors(1:end-1))/2;
m = (Vectors(2:end)+Vectors(1:end-1))/2;
switch(numGaussP)
    case 1
        Gs = 0;
        Gw = 2;
    case 2
        Gs = [-sqrt(1/3), sqrt(1/3)]';
        Gw = [1 ,1 ]';
        m = ones(2,1)*m;
    case 3
        Gs = [-0.774596669241, 0, 0.774596669241]';
        Gw = [5/9, 8/9, 5/9]';
        m = ones(3,1)*m;
    case 4
        Gs = [-0.861136311594053, -0.339981043584856, 0.339981043584856,0.861136311594053]';
        Gw = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454]';
        m = ones(4,1)*m;
    case 5
        Gs = [-0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664]';
        Gw = [0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189]';
        m = ones(5,1)*m;
end
Gs= Gs * h + m;
Gw= Gw * h;
end