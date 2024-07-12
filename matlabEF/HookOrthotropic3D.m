function [D] = HookeOrthotropic3D(E1,E2,E3,nu23,nu31,nu12, ...
                                  G23,G31,G12,mode1,varargin)
% Hooke local operator for orthotropic elasticity
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 30 / 12 / 2004
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 29 / 08 / 2007
%   Ajout base locale
%
% Use Cowin notation for stress and strain symmetric tensors
% same for both : [e11 e22 e33 r2.e23 r2.e31 r2.e12]
%
% Actually it is Cowin notations, not Voigt!

if ~strcmp(mode1,'TRID')
  mode1
  error('mode not correct')
end

% Hooke matrix in local basis
  DI = zeros(3,3);
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

% Global basis if provided
  nin1 = nargin - 10;
  if nin1 == 1
    V1 = varargin{1};
%   Rotation matrix (Voigt notation)
    v1 = V1(:,1);
    v1 = v1 / norm(v1);
    v3 = ProdVect(V1);
    v3 = v3 / norm(v3);
    v2 = ProdVect([v3 v1]);
    v2 = v2 / norm(v2);
    r22 = sqrt(2.)/2.;
%    Q = [StdToVoigt(v1*v1') StdToVoigt(v2*v2') StdToVoigt(v3*v3') ...
%         r22*StdToVoigt(v2*v3' + v3*v2') ...
%         r22*StdToVoigt(v3*v1' + v1*v3') ...
%         r22*StdToVoigt(v1*v2' + v2*v1')];
    Q = [StdToVoigt(v1*v1')'
         StdToVoigt(v2*v2')'
         StdToVoigt(v3*v3')'
         r22*StdToVoigt(v2*v3' + v3*v2')'
         r22*StdToVoigt(v3*v1' + v1*v3')'
         r22*StdToVoigt(v1*v2' + v2*v1')'];
%A VERIFIER    Q' = inv(Q)
%disp('HookOrthotropic3D Should be tested')
    D = inv(Q) * D * Q;
%DD    D = Q' * D * Q;
  elseif (nin1 == 0)
  else
    nin1
    error('Bad number of optional arguments')
  end
