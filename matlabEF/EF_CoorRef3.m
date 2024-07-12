function [idans,Xcor2,shp1,xcor2] = EF_CoorRef3(type1,Xcorel1,xcorel1, ...
                                                xcor1,xcrit1,xcrit2)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 23 / 01 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 12 / 04 / 2003
%   Ajout du TRI6
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 23 / 03 / 2004
%   Modification du critere de proximite
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 02 / 07 / 2006
%   Cas des element 2D dans un espace 3D
% 
% Cherche les coordonnees de reference Xcor2 dans un element donne,
% d'un point donne.
% Au passage, trouve si le point est dans l'element ou non,
% donne la valeur des fonctions de forme en ce point,
% donne le point de l'element le plus proche trouve.
%   
% Entrees
%   type1               : type de l'element geometrique
%   Xcorel1(nbnn,idimr) : coordonnees de reference des noeuds de 
%                         l'element
%   xcorel1(nbnn,idim)  : coordonnees reelles des noeuds de l'element
%   xcor1(1,idim)       : coordonnees reelles du point
%   xcrit1              : critere de proximite absolu (dans espace reel)
%   xcrit2              : critere de proximite relatif (dans espace de ref)
% Sorties
%   idans               : si 0, pas trouve,
%                         si -1, point interieur,
%                         sinon, noeud de l'element
%   Xcor2(1,idimr)      : coordonnees de reference du point le plus proche
%   shp1(nbnn)          : valeurs des fonctions de forme en ce point
%   xcor2(1,idim)       : coordonnees reelles de ce point
%   


[nbnn,idim] = size(xcorel1);
idimr = size(Xcorel1,2);
idans = 0;
Xcor2 = zeros(1,idimr);
shp1  = zeros(1,nbnn);
xcor2 = zeros(1,idim);

% Le point est-il un noeud de l'element ?
% """""""""""""""""""""""""""""""""""""""
% (oui si aucun des ecarts n'est > xcrit1)
isnode1 = find(~any(abs(xcorel1 - repmat(xcor1,nbnn,1)) > xcrit1 , 2));
if (length(isnode1) > 1)
  isnode1
  xcorel1
  error('several nodes of the element are at the same place')
end
if isnode1
  idans = isnode1;
  Xcor2 = Xcorel1(isnode1,:);
  shp1(isnode1) = 1.;
  xcor2 = xcorel1(isnode1,:);
  return;
end


% Si ce n'est pas le cas,
% """""""""""""""""""""""

if ((idim == 3) && (idimr == 2))
%
% Elements 2D plonges dans un espace 3D
% """""""""""""""""""""""""""""""""""""
% On passe en element volumique
  XG1 = sum(xcorel1,1) / nbnn; % Centroid
  IN1 = zeros(idim,idim); % Pseudo-inertial matrix
  for ino1 = 1:nbno
    GM1 = (xcoorg1(ino1,:) - XG1)';
    IN1 = IN1 + ((GM1' * GM1) * eye(idim)) - (GM1 * GM1');
  end
  [N,s,junk] = svd(IN1); % Local basis from Inertial matrix
  I3 = s(1,1); I2 = s(2,2); I1 = s(3,3);
  N3 = N(:,1); N2 = N(:,2); N1 = N(:,3);
  J1 = I2 + I3 - I1;
  J2 = I3 + I1 - I2;
  J3 = I1 + I2 - I3;
% The normal is N3, (N1,N2) = local basis of the plane
  toto = GramS([N3 N2 N1]);
  N3 = toto(:,1); N2 = toto(:,2); N1 = toto(:,3);
% Element 3D applati

  switch type1
    case 'TRI3',
%     on passe en PRI6
      xcorel3D = [xcorel3D - xcrit1*N3' , xcorel3D + xcrit1*N3'];
      type2 = 'PRI6';
      Xcorel3D = EF_CoorRefNod(type2);
      [idans,Xcor2,shp1,xcor2] = EF_CoorRef3D(type2,Xcorel3D,xcorel3D, ...
                                              xcor1,xcrit1,xcrit2);
%     on revient au TRI3
      Xcor2 = Xcor2(1,1:2);
      shp1  = EF_Shape(type1,Xcor2);
      xcor2 = shp1 * xcorel1;

    otherwise,
      type1
      error('element not yet implemented')
  end

elseif ((idim == 3) && (idimr == 3))
%
% Elements massifs 3D
% """""""""""""""""""
  [idans,Xcor2,shp1,xcor2] = EF_CoorRef3D(type1,Xcorel1,xcorel1, ...
                                          xcor1,xcrit1,xcrit2);
elseif ((idim == 2) && (idimr == 2))
%
% Elements massifs 2D
% """""""""""""""""""
  [idans,Xcor2,shp1,xcor2] = EF_CoorRef2D(type1,Xcorel1,xcorel1, ...
                                          xcor1,xcrit1,xcrit2);
end
