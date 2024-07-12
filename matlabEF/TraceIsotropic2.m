function [TR] = TraceIsotropic2(E,nu,mode1)
% Construction de la matrice de prise de la trace des deformations
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 11 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 11 / 2003
%   Modification des composantes 2D ZZ
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 09 / 2007
%   Ajout du mode AXIS
%
% En 2D on utilise la notation de Voith.
%
% En 3D et 2D Deformations Planes, pas besoin de E et nu,
% par contre en 2D Contraintes Planes, il y en a besoin
% (voir HookeIsotropic), mais l'operateur ne fonctionne qu'en lineaire
% en on suppose alors que la composante EPZZ est = 0 ! (laid)
%
% 3D
%                      =   | epsil11 |
%                      =   | epsil22 |
% trace                = TR| epsil33 |
%                      =   | sqrt(2) epsil23 |
%                      =   | sqrt(2) epsil31 |
%                      =   | sqrt(2) epsil12 |

GlobalVar;

% on devrait avoir en 2D (CP et DP) TR = [1. 1. 0. 1.] ??

if strcmp(mode1,liste_mode{1})
% TRID
  TR = [1. 1. 1. 0. 0. 0.];
elseif strcmp(mode1,liste_mode{2})
% COPL Plane stress (Hooke is required)
  D = HookeIsotropic3D(E,nu);
  ddli = [1 2 6 3];
  ddlb = [3];
  TR = [1. 1. 0. 0.] - (D(ddlb,ddlb) \ D(ddlb,ddli));
elseif (strcmp(mode1,'DEPL') || strcmp(mode1,'AXIS'))
% DEPL Plane strain | AXIS Axisymmetric
  TR = [1. 1. 0. 1.];
else
  mode1
  error('Bad mode')
end
