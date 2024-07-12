function [D] = HookePlateS(E,nu,h,mode1)
% Hooke's matrix for a plate model with shear
%
% DUREISSEIX David  LaMCoS                           le 08 / 10 / 2010
%
% Construction de la matrice de Hooke, pour un modele de plaque
% avec cisailement (coief de cisaillement k = 5/6 pour l'instant)
% (en isotrope, sans excentration)
%
% Entrees
%   E		module de Young
%   nu		coefficient de Poisson
%   h		epaisseur
%   mode1	mode d'analyse (DSHE)
% Sorties
%   D		matrice de Hooke
%
% En 3D pour DSHE
% (notation de Cowin : EG12 = sqrt(2) EF12, MG12 = sqrt(2) MF12)
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
% | T1 |      Eh   | G1 |
% | T2 | = k------ | G2 |
%           (1+nu)

H1 = E / (1. - nu^2);
k = 5./6.;
H2 = k*E*h / (1. + nu); % Cowin: no 0.5 coefficient

if strcmp(mode1,'DSHE')
  H = H1 * [1. nu  0.
            nu 1.  0.
            0. 0.  (1 - nu)];
  D = [(h * H)    zeros(3,3)
       zeros(3,3) (h^3 / 12.) * H];
  D = [ D          zeros(6,2)
        zeros(2,6) H2*eye(2,2)];
else
  mode1
  error('Not yet implemented')
end
