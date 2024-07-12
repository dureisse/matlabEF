function [D] = HookeIsotropic3(E,nu,mode1)
% Construction de la matrice de Hooke, en notation de Voith
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 31 / 07 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 01 / 07 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 11 / 2003
%   Modification composantes 2D ZZ
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 09 / 2007
%   Ajout mode AXIS
%
% Entrees
%   E		module de Young
%   nu		coefficient de Poisson
%   mode1	mode d'analyse (TRID,COPL,DEPL,AXIS)
% Sorties
%   D		matrice de Hooke
%
% en 3D (TRID), on retourne D
% | sigma11 |          =   | epsil11 |
% | sigma22 |          =   | epsil22 |
% | sigma33 |          = D | epsil33 |
% | sqrt(2) sigma23 |  =   | sqrt(2) epsil23 |
% | sqrt(2) sigma31 |  =   | sqrt(2) epsil31 |
% | sqrt(2) sigma12 |  =   | sqrt(2) epsil12 |
%
% en 2D (COPL,DEPL,AXIS)
% | sigma11 |          =   | epsil11 |
% | sigma22 |          =   | epsil22 |
% | sqrt(2) sigma12 |  = D | sqrt(2) epsil12 |
% | sigma33 |          =   | epsil33 |
% c'est a dire
% | sigma_plan |         | Dplan  Dhp'  | | epsil_plan |
% | sigma33    |       = | Dhp    Dhors | | epsil33    |
% en deformations planes (DEPL) ou axisymetrique (AXIS),
% on retourne tout le D precedant ;
% en contraintes planes(COPL), on a besoin d'un comportement lineaire,
% (donc de E et nu) pour retourner 
% | sigma_plan |         | Dplan-Dhp'*Dhors^-1*Dhp  0 | | epsil_plan |
% | sigma33    |       = | 0                        0 | | ???????    |

GlobalVar;

D = HookeIsotropic3D(E,nu);

if strcmp(mode1,liste_mode{1})
% TRID

elseif strcmp(mode1,liste_mode{2})
% COPL Plane stress
  ddli = [1 2 6];
  ddlb = [3];
  DD = D(ddli,ddli) - D(ddli,ddlb) * (D(ddlb,ddlb) \ D(ddlb,ddli));
  D = [DD zeros(length(ddli),length(ddlb))
       zeros(length(ddlb),length(ddli)+length(ddlb))];

elseif (strcmp(mode1,'DEPL') || strcmp(mode1,'AXIS'))
% DEPL Plane strain | AXIS Axisymmetric
  ddli = [1 2 6 3];
  ddlb = [4 5];
  D  = D(ddli,ddli);

else
  mode1
  error('Bad mode')
end
