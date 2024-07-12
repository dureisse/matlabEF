function [D] = HookeAnisotropic3D(YG1,YG2,YG3,NU12,NU23,NU13,G12,G23,G13, ...
                                  V1X,V1Y,V1Z,V2X,V2Y,V2Z)
% Hooke local operator for anisotropic elasticity
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 29 / 08 / 2007
%
% Use Voigt notation for stress and strain symmetric tensors
% same for both : [e11 e22 e33 r2.e23 r2.e31 r2.e12]

% Hooke matrix in local basis
error('to be done')
  DI = zeros(6,6);

  DI(1,1) = 1. ./ E1;
  DI(1,2) = -nu12 ./ E1;
  DI(1,3) = -nu31 ./ E3;
  DI(2,1) = -nu12 ./ E1;
  DI(2,2) = 1. ./ E2;
  DI(2,3) = -nu23 ./ E2;
  DI(3,1) = -nu31 ./ E3;
  DI(3,2) = -nu23 ./ E2;
  DI(3,3) = 1. ./ E3;

  D = zeros(6,6);
  D(1:3,1:3) = inv(DI); clear DI;
  D(4,4) = 2.*G23;
  D(5,5) = 2.*G31;
  D(6,6) = 2.*G12;
