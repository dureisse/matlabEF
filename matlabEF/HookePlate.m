function [D] = HookePlate(E,nu,h,mode1)
% Hooke's matrix for a plate model
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 07 / 2004
%
% Construction de la matrice de Hooke, pour un modele de plaque
% (Kirchhoff en dimension 3 pour l'instant...)
% (en isotrope, sans excentration)
%
% Entrees
%   E		module de Young
%   nu		coefficient de Poisson
%   h		epaisseur
%   mode1	mode d'analyse (DKIR)
% Sorties
%   D		matrice de Hooke
%
% En 3D pour DKIR
% (notation de Voith : EG12 = sqrt(2) EF12, MG12 = sqrt(2) MF12)
% | EF11 |         |         0 0 0   | | EP11 |
% | EF22 |      E  |  h H    0 0 0   | | EP22 |
% | EG12 | = ------|         0 0 0   |*| GA12 |
% | MF11 |   1-nu^2| 0 0 0           | | CP11 |
% | MF22 |         | 0 0 0  h^3/12 H | | CP22 |
% | MG12 |         | 0 0 0           | | CG12 |
%
% avec     | 1  nu 0        |
%      H = | nu 1  0        |
%          | 0  0  (1-nu)/4 |
%

H1 = E / (1. - nu^2);

if strcmp(mode1,'DKIR')
  H = H1 * [1. nu  0.
            nu 1.  0.
            0. 0.  (1 - nu)];
  D = [(h * H)    zeros(3,3)
       zeros(3,3) (h^3 / 12.) * H];
else
  mode1
  error('Not yet implemented')
end
