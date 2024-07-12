function [idans,Xcor2,shp1,xcor2] = EF_CoorRef2D(type1,Xcorel1,xcorel1, ...
                                                xcor1,xcrit1,xcrit2)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 23 / 01 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 12 / 04 / 2003
%   Ajout du TRI6
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 23 / 03 / 2004
%   Modification du critere de proximite
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

% On sous-decoupe l'element en triangles, et on cherche les coordonnees
% barycentriques dans les triangles

% Barycentre de reference (coordonnees de ref Xp0 et reelles xp0)
Xp0 = sum(Xcorel1,1) / nbnn;
shp1 = EF_Shape(type1,Xp0);
xp0 = shp1 * xcorel1;

if strcmp(type1,'TRI3')

  Xp0 = Xcorel1(1,:);
  xp0 = xcorel1(1,:);
  clear ial; ial(1) = 2; ial(2) = 3;
  [idans,al] = Efbary(xp0,xcorel1',ial,xcor1,xcrit2);

elseif strcmp(type1,'QUA4')

  if ~idans
    clear ial; ial(1) = 1; ial(2) = 2;
    [idans,al] = Efbary(xp0,xcorel1',ial,xcor1,xcrit2);
  end
  if ~idans
    clear ial; ial(1) = 2; ial(2) = 3;
    [idans,al] = Efbary(xp0,xcorel1',ial,xcor1,xcrit2);
  end
  if ~idans
    clear ial; ial(1) = 3; ial(2) = 4;
    [idans,al] = Efbary(xp0,xcorel1',ial,xcor1,xcrit2);
  end
  if ~idans
    clear ial; ial(1) = 4; ial(2) = 1;
    [idans,al] = Efbary(xp0,xcorel1',ial,xcor1,xcrit2);
  end

elseif strcmp(type1,'TRI6')

  if ~idans
    clear ial; ial(1) = 1; ial(2) = 4;
    [idans,al] = Efbary(xp0,xcorel1',ial,xcor1,xcrit2);
  end
  if ~idans
    clear ial; ial(1) = 4; ial(2) = 2;
    [idans,al] = Efbary(xp0,xcorel1',ial,xcor1,xcrit2);
  end
  if ~idans
    clear ial; ial(1) = 2; ial(2) = 5;
    [idans,al] = Efbary(xp0,xcorel1',ial,xcor1,xcrit2);
  end
  if ~idans
    clear ial; ial(1) = 5; ial(2) = 3;
    [idans,al] = Efbary(xp0,xcorel1',ial,xcor1,xcrit2);
  end
  if ~idans
    clear ial; ial(1) = 3; ial(2) = 6;
    [idans,al] = Efbary(xp0,xcorel1',ial,xcor1,xcrit2);
  end
  if ~idans
    clear ial; ial(1) = 6; ial(2) = 1;
    [idans,al] = Efbary(xp0,xcorel1',ial,xcor1,xcrit2);
  end

else
  type1
  error('element not implemented yet')
end


if idans
  idans = -1;
% Rem : avec min(ial), max(ial) on peut savoir si on est en un noeud,
% sur le bord, sur la face !

% Le point est a l'interieur : on cherche les coordonnees de reference
% """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
% on initialise la position de reference
  Xcor2 = al(1) * Xp0;
  for i = 1:idim
    Xcor2 = Xcor2 + al(i+1) * Xcorel1(ial(i),:);
  end
%%  [Xcor2,dist1] = fminsearch(@Efdist,Xcor2,[],xcor1,type1,xcorel1);
  options = optimset('Display','notify', ...
                     'MaxFunEvals','200*numberofvariables', ...
		     'MaxIter','200*numberofvariables', ...
		     'TolFun',min(1.e-4,xcrit1), ...
		     'TolX',min(1.e-4,xcrit2));
  [Xcor2,dist1] = fminsearch('Efdist',Xcor2,options,xcor1,type1,xcorel1);
if (dist1 > xcrit1)
   dist1
  xcrit1
  disp('Point mal trouve')
  keyboard
  error('Point mal trouve')
end
  shp1 = EF_Shape(type1,Xcor2);
  xcor2 = shp1 * xcorel1;

end





% =====================================================================
function [idans,AL] = Efbary(XPP,XCOR3,IAL,XCOR2,xcrit2)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 23 / 01 / 2003
%
% Recherche des coordonnees barycentriques reelles
%
% Entrees
%   XPP(idim)        : coordonnees du premier point
%   XCOR3(idim,nbnn) : coordonnees reelles des points
%   IAL(idim)        : ordre des noeuds qui interviennent
%   XCOR2(idim)      : coordonnees du point cherche
%   xcrit2           : critere de proximite relatif
% Sorties
%   idans            : si 1, a l'interieur ; si 0, a l'exterieur
%   AL(idim+1)       : coordonnees barycentriques
%
% Remarque : le critere de proximite est au sens large

idim = size(XCOR3,1);
NBNG = idim + 1;
AL = zeros(1,NBNG);
if (idim == 2)
             X1=XPP(1);
             X2=XCOR3(1,IAL(1));
             X3=XCOR3(1,IAL(2));
             Y1=XPP(2);
             Y2=XCOR3(2,IAL(1));
             Y3=XCOR3(2,IAL(2));
             X=XCOR2(1);
             Y=XCOR2(2);
% -
             DETAM=X1*Y2+X2*Y3+X3*Y1-Y1*X2-Y2*X3-Y3*X1;
             A(1)=X2*Y3-X3*Y2;
             A(2)=X3*Y1-X1*Y3;
             A(3)=X1*Y2-X2*Y1;
             B(1)=Y2-Y3;
             B(2)=Y3-Y1;
             B(3)=Y1-Y2;
             C(1)=X3-X2;
             C(2)=X1-X3;
             C(3)=X2-X1;
             for IK=1:NBNG
               AL(IK)=(A(IK)+B(IK)*X+C(IK)*Y)/DETAM;
             end
             AL(4)=1.;
elseif (idim == 3)
             X1=XPP(1);
             X2=XCOR3(1,IAL(1));
             X3=XCOR3(1,IAL(2));
             X4=XCOR3(1,IAL(3));
             Y1=XPP(2);
             Y2=XCOR3(2,IAL(1));
             Y3=XCOR3(2,IAL(2));
             Y4=XCOR3(2,IAL(3));
             Z1=XPP(3);
             Z2=XCOR3(3,IAL(1));
             Z3=XCOR3(3,IAL(2));
             Z4=XCOR3(3,IAL(3));
             X=XCOR2(1);
             Y=XCOR2(2);
             Z=XCOR2(3);
% -
  DETAM=X2*Y3*Z4+X3*Y4*Z2+X4*Y2*Z3-X4*Y3*Z2-X2*Y4*Z3-X3*Y2*Z4-X1*Y3* ...
  Z4-X3*Y4*Z1-X4*Y1*Z3+X4*Y3*Z1+X3*Y1*Z4+X1*Y4*Z3+X1*Y2*Z4+X4*Y1*Z2+ ...
  X2*Y4*Z1-X4*Y2*Z1-X2*Y1*Z4-X1*Y4*Z2-X1*Y2*Z3-X3*Y1*Z2-X2*Y3*Z1+X3* ...
  Y2*Z1+X2*Y1*Z3+X1*Y3*Z2;
             A(1)=X2*Y3*Z4+X3*Y4*Z2+X4*Y2*Z3-X4*Y3*Z2-X2*Y4*Z3-X3*Y2*Z4;
             A(2)=X4*Y3*Z1+X3*Y1*Z4+X1*Y4*Z3-X1*Y3*Z4-X3*Y4*Z1-X4*Y1*Z3;
             A(3)=X1*Y2*Z4+X4*Y1*Z2+X2*Y4*Z1-X4*Y2*Z1-X2*Y1*Z4-X1*Y4*Z2;
             A(4)=X3*Y2*Z1+X2*Y1*Z3+X1*Y3*Z2-X1*Y2*Z3-X3*Y1*Z2-X2*Y3*Z1;
             B(1)=Y4*Z3-Y3*Z4+Y2*Z4-Y4*Z2+Y3*Z2-Y2*Z3;
             B(2)=Y3*Z4-Y4*Z3+Y4*Z1-Y1*Z4+Y1*Z3-Y3*Z1;
             B(3)=Y4*Z2-Y2*Z4+Y1*Z4-Y4*Z1+Y2*Z1-Y1*Z2;
             B(4)=Y2*Z3-Y3*Z2+Y3*Z1-Y1*Z3+Y1*Z2-Y2*Z1;
             C(1)=X3*Z4-X4*Z3+X4*Z2-X2*Z4+X2*Z3-X3*Z2;
             C(2)=X4*Z3-X3*Z4+X1*Z4-X4*Z1+X3*Z1-X1*Z3;
             C(3)=X2*Z4-X4*Z2+X4*Z1-X1*Z4+X1*Z2-X2*Z1;
             C(4)=X3*Z2-X2*Z3+X1*Z3-X3*Z1+X2*Z1-X1*Z2;
             D(1)=X4*Y3-X3*Y4+X2*Y4-X4*Y2+X3*Y2-X2*Y3;
             D(2)=X3*Y4-X4*Y3+X4*Y1-X1*Y4+X1*Y3-X3*Y1;
             D(3)=X4*Y2-X2*Y4+X1*Y4-X4*Y1+X2*Y1-X1*Y2;
             D(4)=X2*Y3-X3*Y2+X3*Y1-X1*Y3+X1*Y2-X2*Y1;
             for II=1:NBNG
               AL(II)=(A(II)+B(II)*X+C(II)*Y+D(II)*Z)/DETAM;
             end 
else
  idim
  error('bad idim')
end
if all(AL > -1.*xcrit2)
% Le point est interne a l'element
  AL = (abs(AL) > xcrit2) .* AL;
  idans = 1;
else
%AL
%-1.*xcrit2
  idans = 0;
end
