function [D] = HookeIsotropic3D(E,nu)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 11 / 2002

% Hooke local operator for isotropic elasticity

  D = zeros(6,6);
  DMU = E ./ (1. + nu);
  LAM = (DMU .* nu) ./ (1. - 2. * nu);
  LPD = LAM + DMU;
  D(1,1) = LPD;
  D(1,2) = LAM;
  D(1,3) = LAM;
  D(2,1) = LAM;
  D(2,2) = LPD;
  D(2,3) = LAM;
  D(3,1) = LAM;
  D(3,2) = LAM;
  D(3,3) = LPD;
  D(4,4) = DMU;
  D(5,5) = DMU;
  D(6,6) = DMU;
