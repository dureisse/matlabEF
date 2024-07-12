function [dist1] = Efdist(Xcor1,xcor2,type1,xcorel1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 23 / 01 / 2003
%
% Calcule la distance reelle entre le point de coordonnees reelles
% xcor2, et le point dont les coordonnees de references dans 
% l'element fini sont Xcor1
%
% Entrees
%  Xcor1(idimr)       : coordonnees de reference du point teste
%  xcor2(idim)        : coordonnees reelles du point cherche
%  type1              : type de l'element geometrique
%  xcorel1(nbnn,idim) : coordonnees reelles des noeuds de l'element
% Sortie
%  dist1              : ecart entre le point cherche et le point teste
%                       (en coordonnees reelles)

% Coordonnees reelles du point teste xcor1
  shp1 = EF_Shape(type1,Xcor1);
  xcor1 = shp1 * xcorel1;

% Distance
  dist1 = norm(xcor1 - xcor2);
