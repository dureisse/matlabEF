function [val1] = IntegralTheta(t1,f1,theta1,varargin)
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 27 / 05 / 2003
%
% Calcul de l'integrale d'une fonction f1 sur la grille en temps t1
% par la theta-methode
% (voir IntegrateTheta)
%
% Entrees
%   t1(ntmps)		Liste des piquets de temps
%   f1(npts,ntmps)	Liste des valeurs de la fonction a integrer
%   theta1		Parametre de la methode d'integration en temps
% Entrees optionnelles
%   weight(ntmps)	Liste des ponderations en temps eventuelle
% Sorties
%   val1(npts)		Valeur de l'integrale en temps

[npts,ntmps] = size(f1);

nin1 = nargin-3;
if (nin1 == 0)
  weight = ones(1,ntmps);
elseif (nin1 == 1)
  weight = varargin{1};
else
  nin1
  error('Bad number of arguments')
end

it1 = 1;
  val1 = zeros(npts,1) * weight(it1);

for it1 = 2:ntmps
  hi = t1(it1) - t1(it1-1);
  val1 = val1 + hi*(theta1*f1(:,it1)*weight(it1) + ...
                    (1.-theta1)*f1(:,it1-1)*weight(it1-1));
end
