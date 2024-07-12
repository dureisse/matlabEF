function [D] = HookeBeam(E,nu,S,IZ,mode1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 09 / 12 / 2003
%
% Construction de la matrice de Hooke, pour un modele de poutre
% (Timoshenko ou Euler-Bernoulli, en dimension 2 pour l'instant...)
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
% | EFF1 |   | ES 0  0    | | EPS1 |
% | EFF2 | = | 0  GS 0    |*| GA12 |
% | MOM3 | = | 0  0  E.IZ | | C3   |
% En 2D pour POUT
% | EFF1 |   | ES 0    | | EPS1 |
% | MOM3 | = | 0  E.IZ | | C3   |

G = 0.5 * E / (1. + nu);
if strcmp(mode1,'TIMO')
  D = [E*S 0.  0.
       0.  G*S 0.
       0.  0.  E*IZ];
elseif strcmp(mode1,'POUT')
  D = [E*S 0.
       0.  E*IZ];
else
  mode1
  error('Not yet implemented')
end
