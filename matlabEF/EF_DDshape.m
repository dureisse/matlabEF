function [ddshp1] = EF_DDshape(type1,X);
% Second derivative of finite element shape functions
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 11 / 07 / 2004
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 02 / 07 / 2006
%  Ajout TET4, PRI6
% DUREISSEIX David  LaMCoS                          le 08 / 10 / 2010
%  Ajout DSQ, DST
% 
% Retourne les valeurs des derivees secondes (par rapport aux
% coordonnees de reference) des fonctions de forme d'un element
% geometrique de type type1 au point dont les coordonnees de
% reference sont X.
%   
% Entrees
%   type1               : type de l'element geometrique
%   X(idimr)            : coordonnees de reference du point
% Sorties
%   ddshp1(idimr,idimr,nbnn)   : valeurs des derivees secondes en ce point

switch type1
  case 'SEG2',
    ddshp1 = zeros(1,1,2);

  case 'TRI3',
    ddshp1 = zeros(2,2,3);

  case 'QUA4',
    ddshp1 = zeros(2,2,4);
    ddshp1(:,:,1) = 0.25 * [0.  1. ;  1. 0.];
    ddshp1(:,:,2) = 0.25 * [0. -1. ; -1. 0.];
    ddshp1(:,:,3) = 0.25 * [0.  1. ;  1. 0.];
    ddshp1(:,:,4) = 0.25 * [0. -1. ; -1. 0.];

  case {'DKQ','DSQ'}
    ddshp1 = zeros(2,2,4+4);
    ddshp1(:,:,1) = 0.25 * [0.  1. ;  1. 0.];
    ddshp1(:,:,2) = 0.25 * [0. -1. ; -1. 0.];
    ddshp1(:,:,3) = 0.25 * [0.  1. ;  1. 0.];
    ddshp1(:,:,4) = 0.25 * [0. -1. ; -1. 0.];
    ddshp1(:,:,5) = [(X(2)-1.) X(1) ; X(1) 0.];
    ddshp1(:,:,6) = [0. -X(2) ; -X(2) -1.-X(1)];
    ddshp1(:,:,7) = [-1.-X(2) -X(1) ; -X(1) 0.];
    ddshp1(:,:,8) = [0. X(2) ; X(2) (X(1)-1.)];

  case {'DKT','DST'}
    ddshp1 = zeros(2,2,3+3);
    ddshp1(:,:,1:3) = zeros(2,2,3);
    ddshp1(:,:,4) = [-8. -4. ; -4. 0.];
    ddshp1(:,:,5) = [0. 4. ; 4. 0.];
    ddshp1(:,:,6) = [0. -4. ; -4. -8.];

  case 'TET4',
    ddshp1 = zeros(3,3,4);

  case 'PRI6',
    ddshp1 = zeros(3,3,6);
    ddshp1(:,:,1) = 0.5*[0. 0.  1. ; 0. 0.  1. ;  1.  1. 0.];
    ddshp1(:,:,2) = 0.5*[0. 0. -1. ; 0. 0.  0. ; -1.  0. 0.];
    ddshp1(:,:,3) = 0.5*[0. 0.  0. ; 0. 0. -1. ;  0. -1. 0.];
    ddshp1(:,:,4) = 0.5*[0. 0. -1. ; 0. 0. -1. ; -1. -1. 0.];
    ddshp1(:,:,5) = 0.5*[0. 0.  1. ; 0. 0.  0. ;  1.  0. 0.];
    ddshp1(:,:,6) = 0.5*[0. 0.  0. ; 0. 0.  1. ;  0.  1. 0.];

  otherwise,
    type1
    error('Type of element not implemented yet')
end

