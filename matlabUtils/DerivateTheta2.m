function [Df1] = DerivateTheta2(t1,f1,theta1,xval0)
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 12 / 12 / 2002
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 25 / 05 / 2003
%   Possibilite d'avoir un champ f1(:,:, ... ,t)
%
% Df1 = DerivateTheta2(t1,f1,theta1,xval0)
% Derivation d'une fonction f1 sur la grille en temps t1
% par la theta-methode
% Le resultat est la fonction Df1 sur la grille en temps t1
% On a besoin de la valeur initiale de la derivee Df1(t=0) = xval0

Df1 = 0. * f1;

str = repmat(':,',1,ndims(f1)-1);

t = 1;
% Df1(::,1) = xval0;
eval(['Df1(' str int2str(t) ') = xval0;'])

for t = 2:length(t1)
  hi = t1(t) - t1(t-1);
%  Df1(::,it1) = (f1(::,it1) - f1(::,it1-1))/(theta1*hi) - ...
%             ((1.-theta1)/theta1)*Df1(::,it1-1);
eval(['Df1(' str int2str(t) ') = (f1(' str int2str(t) ') - f1(' str int2str(t-1) '))/(theta1*hi) - ((1.-theta1)/theta1)*Df1(' str int2str(t-1) ');'])
end
