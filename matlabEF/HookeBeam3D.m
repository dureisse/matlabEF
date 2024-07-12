function [D] = HookeBeam3D(E,nu,S,I,I0,mode1)
% Hooke matrix for a 3D beam model (notice that I2=I3 up to now!)
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 17 / 01 / 2006
%
% Construction de la matrice de Hooke, pour un modele de poutre
% (Timoshenko ou Euler-Bernoulli, en dimension 3)
%
% Entrees
%   E		module de Young
%   nu		coefficient de Poisson
%   S		section
%   IZ		moment de surface de la section droite
%   mode1	mode d'analyse (TIMO,POUT)
% Sorties
%   D		matrice de Hooke
%
% En 2D pour TIMO
% Not yet done for 3D
% | EFF1 |   | ES 0  0    | | EPS1 |
% | EFF2 | = | 0  GS 0    |*| GA12 |
% | MOM3 | = | 0  0  E.IZ | | C3   |
% En 3D pour POUT (MOM1 correpond a la torsion)
% | EFF1 |   | ES 0    0   0   | | EPS1 |
% | MOM1 | = | 0  G.I0 0   0   | | C1   |
% | MOM2 | = | 0  0    E.I 0   | | C2   |
% | MOM3 | = | 0  0    0   E.I | | C3   |

G = 0.5 * E / (1. + nu);
if strcmp(mode1,'TIMO')
  error('to be done for TIMO')
  D = [E*S 0.  0.
       0.  G*S 0.
       0.  0.  E*IZ];
elseif strcmp(mode1,'POUT')
  D = [E*S 0.   0.  0.
       0.  G*I0 0.  0.
       0.  0.   E*I 0.
       0.  0.   0.  E*I];
else
  mode1
  error('Not yet implemented')
end
