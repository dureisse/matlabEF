function [N,I,G] = GeometricForm(xcoor,varargin)
% Try to identify the form of a cloud of nodes
%
% DUREISSEIX David  L.M.G.C.  SYSTEMES MULTICONTACTS le 15 / 12 / 2006
%
% Inputs
%   xcoor(nbno,3)	coordinates of the nodes
% Optional inputs
%   w(nbno)		weights (default: 1)
% Ouputs
%   N(3,3)		local basis N1=N(:,1) N2=N(:,2) N3=N(:,3)
%   I(3)		associated proper values (dimensions)
%   G(3)		centroid coordinates

  nin1 = nargin - 1;
  switch nin1
    case 0,
      w = 1;
    case 1,
      w = varargin{1};
      error('Weight not yet done')
    otherwise,
      nin1
      error('Bad number of optional input arguments')
  end

  idim = size(xcoor,2);
  if (idim ~= 3)
    idim
    error('Only in 3D')
  end

  nbno = size(xcoor,1);

% Centroid
  G = sum(xcoor,1) / nbno;

% Pseudo-inertial matrix
  IN1 = zeros(idim,idim);
  for ino1 = 1:nbno
    GM1 = (xcoor(ino1,:) - G)';
    IN1 = IN1 + ((GM1' * GM1) * eye(idim)) - (GM1 * GM1');
  end

% Local basis from Inertial matrix
  [N,s,junk] = svd(IN1);
  I3 = s(1,1); I2 = s(2,2); I1 = s(3,3);
  N3 = N(:,1); N2 = N(:,2); N1 = N(:,3);
  J1 = I2 + I3 - I1;
  J2 = I3 + I1 - I2;
  J3 = I1 + I2 - I3;

% In the case of a plane, the normal is N3, (N1,N2) = local basis of the plane
  toto = GramS([N3 N2 N1]);
  N3 = toto(:,1); N2 = toto(:,2); N1 = toto(:,3);

  N = [N1 N2 N3];
  I = [J1 J2 J3];
