function [Mjaco,Jaco] = LocalJaco2(dphix,xcoor)

% Matrice Jacobienne et Jacobien de la transformation
Mjaco = dphix * xcoor;
[idimr,idim] = size(Mjaco);
if (idim == idimr)
  Jaco = det(Mjaco);
elseif ((idim == 2) && (idimr == 1))
% Element 1D dans un espace 2D
  Jaco = norm(Mjaco);
elseif ((idim == 3) && (idimr == 1))
% Element 1D dans un espace 3D
  Jaco = norm(Mjaco);
elseif ((idim == 3) && (idimr == 2))
% Element 2D dans un espace 3D
  Jaco = norm(ProdVect(Mjaco'));
else
  idim
  idimr
  error('plongement pas prevu')
end
