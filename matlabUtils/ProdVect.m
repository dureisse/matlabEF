function [resu] = ProdVect(vects)
%
% Produit vectoriel resu(n,1) des vecteurs de vects(n,nvect)
% Il faut nvect = n-1 pour avoir un vecteur en sortie, 
% si nvect = n, on a la projection du vecteur resultat sur la normale
% en 2D (determinant), le produit mixte en 3D

[n,nvect] = size(vects);
%if ~(nvect == n-1)
%  nvect
%  n
%  error('bad dim')
%end

if (nvect == n-1)
  resu = zeros(n,1);
  if (n == 2)
    resu(1,1) = -vects(2,1);
    resu(2,1) = vects(1,1);
  elseif (n == 3)
    resu(1,1) = vects(2,1)*vects(3,2)-vects(3,1)*vects(2,2);
    resu(2,1) = vects(3,1)*vects(1,2)-vects(1,1)*vects(3,2);
    resu(3,1) = vects(1,1)*vects(2,2)-vects(2,1)*vects(1,2);
  else
    n
    error('dimension pas prevue')
  end
elseif (nvect == n)
  resu = det(vects);
else
  error('nombre de vecteurs pas prevu')
end
