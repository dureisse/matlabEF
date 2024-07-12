function [Dg] = HookeLocalToGlobal(Dl,Q);
% Rotation for Hooke local matrix
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 30 / 12 / 2004
%
% Input
%   Dl(n,n)		Hooke matrix in local coordinates
%   Q(idim,idim)	Rotation matrix (coordinates of local vectors
%                       in global basis)
% Output
%   Dg(n,n)		Hooke matrix in global coordinates
%
% Use Voigt notation for stress and strain symmetric tensors
% same for both : [e11 e22 e33 r2.e23 r2.e31 r2.e12] for 3D
%                 [e11 e22 r2.e12 e33] for 2D

  idim = size(Q,1);
  n = size(Dl,1);

  r2 = sqrt(2.);

% Rotation for Voigt notation
% """""""""""""""""""""""""""
  if (idim == 2)
%   2D comme cas particulier de 3D
    Qnew = eye(6,6);
    Qnew([1 2 6 3],[1 2 6 3]) = Q;
    Q = Qnew;
  end

% 3D case
  QQ = zeros(6,6);
  for i = 1:3
    for j = 1:3
      QQ(i,j) = Q(i,j)^2;
    end
    QQ(i,4) = r2 * Q(i,2) * Q(i,3);
    QQ(i,5) = r2 * Q(i,3) * Q(i,1);
    QQ(i,6) = r2 * Q(i,1) * Q(i,2);

    QQ(4,i) = r2 * Q(2,i) * Q(3,i);
    QQ(5,i) = r2 * Q(3,i) * Q(1,i);
    QQ(6,i) = r2 * Q(1,i) * Q(2,i);
  end

  QQ(4,4) = Q(2,2) * Q(3,3) + Q(2,3) * Q(3,2);
  QQ(4,5) = Q(2,3) * Q(3,1) + Q(2,1) * Q(3,3);
  QQ(4,6) = Q(2,1) * Q(3,2) + Q(2,2) * Q(3,1);

  QQ(5,4) = Q(3,2) * Q(1,3) + Q(3,3) * Q(1,2);
  QQ(5,5) = Q(3,3) * Q(1,1) + Q(3,1) * Q(1,3);
  QQ(5,6) = Q(3,1) * Q(1,2) + Q(3,2) * Q(1,1);

  QQ(6,4) = Q(1,2) * Q(2,3) + Q(1,3) * Q(2,2);
  QQ(6,5) = Q(1,3) * Q(2,1) + Q(1,1) * Q(2,3);
  QQ(6,6) = Q(1,1) * Q(2,2) + Q(1,2) * Q(2,1);

% Back to 2D
  if (idim == 2)
    QQold = QQ([1 2 6 3],[1 2 6 3]);
    QQ = QQold;
  end

% Applying rotation
% """""""""""""""""
  Dg = QQ * Dl * QQ';

  clear QQ QQold Qnew;
