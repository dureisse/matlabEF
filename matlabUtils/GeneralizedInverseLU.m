function [Ap,R,A] = GeneralizedInverseLU(A,xcrit1);
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 08 / 05 / 00
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 03 / 01 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 13 / 08 / 2006
%
% Inversion d'une matrice carree non symetrique A par la methode LU,
% sans pivotage, mais avec recuperation des modes rigides.
% A contient en sortie sa factorisation LU generalisee
% (si mode rigide detecte en position I, D(I) est mis a 1)
% 
% Input
%   A(N,N)      Matrice a inverser
%   xcrit1	Precision pour la detection des modes rigides
% Output
%   Ap(N,N)     Inverse generalisee
%   R(N,NRIG)   Modes rigides




% Factorisation LU
[L,U] = lu(A); % L*U = A

% Search zeros on diagonal of U
ds = abs(diag(U));
s0 = max(ds);
ds = full(ds) / s0;
[ss,ind] = sort(ds); % ss = ds(ind)
k = 0; while ss(k+1) < xcrit1
  k = k + 1;
end
lpivot = ind(1:k)';
%disp([' Kernel size ' int2str(k)])
nddl = size(A,2);
lreg = setdiff([1:nddl],lpivot);

% Fill in R
R = zeros(nddl,k);
Arr = A(lreg,lreg);
Ari = A(lreg,lpivot);
R(lpivot,:) = eye(k);
R(lreg,:) = - Arr \ (Ari * eye(k));

% Generalised Inverse
Ap = pinv(A,xcrit1);
%Pr = R * ((R'*R) \ R');
%Ap = Ap - Pr * Ap * Pr;
