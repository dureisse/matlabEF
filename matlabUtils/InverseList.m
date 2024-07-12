function [liste_inv] = InverseList(liste,nmax)
% Inverse une liste d'indirection
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 01 / 08 / 2002
% 
% liste_inv = InverseList(liste,nmax);
% Input 
%   liste(nb) : liste d'indirection
%   nmax : valeur maximale des numeros
% Output
%   liste(nmax) : liste inverse
%   si liste(i)=j (0< j <= nmax), alors liste_inv(j)=i
%   si j n'est pas dans liste, alors liste_inv(j)=0
%   si liste(i)=0, i ne sera pas dans la liste_inv

liste_inv = zeros(1,nmax);
for i = 1:size(liste,2)
  if liste(i)
    liste_inv(liste(i)) = i;
  end
end
