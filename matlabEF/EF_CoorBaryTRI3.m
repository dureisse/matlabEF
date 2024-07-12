function [T1,kerr] = EF_CoorBaryTRI3(xpoly1,xcoor2,xprec);
% Find barycentric coordinates of nodes within a TRI3 element
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 17 / 03 / 2003
%
% Inputs
%   xpoly1(3,2)		Coordinates of the nodes of the TRI3 element
%   xcoor2(nbno2,2)	Coordinates of the target nodes
%   xprec		Precision (cond)
% Outputs
%   T1(nbno2,3)		Projection operator = barycentric coordinates
%   kerr		Error indicator

kerr = 0;
nbno = size(xpoly1,1);
idim = size(xpoly1,2);
nbno2 = size(xcoor2,1);

if (nbno ~= 3) | (idim ~= 2) | (size(xcoor2,2) ~= idim)
  nbno
  idim
  size(xcoor2,2)
  kerr = 1
  error('Erroneous inputs')
end

% Reference coordinates beta
N = [-1 -1 ; 1 0 ; 0 1];
M = N'*xpoly1*xpoly1'*N;
% test for bad scaling
s=svd(M);
if (s(2)/s(1)) < xprec
  disp('EF_CoorBaryTRI3: Warning: Bad conditionning')
  T1 = zeros(nbno2,3);
  kerr = 1;
  return
end
u = xcoor2 - repmat(xpoly1(1,:),nbno2,1);
beta = M \ (N'*xpoly1*u');

% Barycentric coordinates
alpha1 = 1 - sum(beta,1);
T1 = [alpha1' beta'];
