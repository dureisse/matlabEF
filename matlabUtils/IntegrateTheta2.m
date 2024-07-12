function [f1] = IntegrateTheta2(t1,Df1,theta1,xval0)
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 12 / 12 / 2002
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 10 / 06 / 2003
%   Possibilite d'avoir un champ Df1(:,:, ... ,t)
%
% Integration d'une fonction Df1 sur la grille en temps t1
% par la theta-methode
% Le resultat est la fonction f1 sur la grille en temps t1
% On a besoin de la valeur initiale f1(t=0) = xval0


f1 = 0. * Df1;
str = repmat(':,',1,ndims(Df1)-1);

t = 1;
% f1(::,t) = xval0;
  eval(['f1(' str int2str(t) ') = xval0;'])

for t = 2:length(t1)
  hi = t1(t) - t1(t-1);
%  f1(::,t) = f1(::,t-1) + hi*(theta1*Df1(::,t) + ...
%               (1.-theta1)*Df1(::,t-1));
  eval(['f1(' str int2str(t) ') = f1(' str int2str(t-1) ...
        ') + hi*(theta1*Df1(' ...
        str int2str(t) ') + (1.-theta1)*Df1(' str int2str(t-1) '));'])
end
