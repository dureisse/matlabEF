function [Xcorel1] = EF_CoorRefNod(type1);
% Give natural coordinates Xcorel1 of nodes for element type type1
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 23 / 01 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 09 / 04 / 2003
%   Ajout TRI6
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 21 / 07 / 2003
%   Ajout SEG2
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 15 / 05 / 2005
%   Ajout RAC2
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 16 / 09 / 2005
%   Ajout RAC3
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 02 / 07 / 2006
%   Ajout PRI6, TET4
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 04 / 08 / 2006
%   Ajout CUB8
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 16 / 12 / 2006
%   Ajout POI1
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 03 / 03 / 2007
%   Ajout QUA8, SEG3, PR15
%
% Retourne les coordonnees de reference Xcorel1 des noeuds d'un
% element geometrique de type type1
%   
% Entrees
%   type1               : type de l'element geometrique
% Sorties
%   Xcorel1(nbnn,idimr) : coordonnees de reference des noeuds
%
%     SEG2       RAC2      RAC3
%     1-----2    1----2    1--2--3
%    -1  0  1    0    1   -1  0  1
%
%     SEG3
%     1--3--2
%    -1  0  1
%
%     TRI3
%   1 3
%     |\
%     | \
%     1--2
%   0    1
%
%     QUA4
%   1 3----4
%     |    |
%     |    |
%  -1 1----2
%    -1    1
%
%     QUA8
%   1 3--7--4
%     |     |
%     8     6
%     |     |
%  -1 1--5--2
%    -1    1
%
%     TRI6
%   1 3
%     |\
%     6 5
%     |  \
%     1-4-2
%   0     1
%
%     PRI6
%         0  1
%         4--6
%        /|  |
%     1 5 |  |
%       | 1--3
%       |/
%    -1 2
%
%     TET4
%       1 4
%         |
%         |
%         1--3
%        /   1
%     1 2
%
%     TE10
%       comme gmsh
%
%     CUB8

switch type1
  case 'POI1'
    Xcorel1 = [0.];  % par convention
  case 'SEG2'
    Xcorel1 = [-1.
                1.];
  case 'SEG3'
    Xcorel1 = [-1.
                1.
                0.];
  case 'TRI3',
    Xcorel1 = [0. 0.
               1. 0.
               0. 1.];
  case 'QUA4',
    Xcorel1 = [-1. -1.
                1. -1.
                1. 1.
               -1. 1.];
  case 'QUA8',
    Xcorel1 = [-1. -1.
                1. -1.
                1. 1.
               -1. 1.
                0. -1.
                1. 0.
                0. 1.
               -1. 0.];
  case 'TRI6',
    Xcorel1 = [0.  0.
               1.  0.
               0.  1.
               0.5 0.
               0.5 0.5
               0.  0.5];
  case 'RAC2',
%   Only 2 nodes of 4 taken into account
    Xcorel1 = [0.
               1.];
  case 'RAC3',
%   Only 3 nodes of 6 taken into account
    Xcorel1 = [-1.
                0.
                1.];
  case 'PRI6',
    Xcorel1 = [0. 0. -1.
               1. 0. -1.
               0. 1. -1.
               0. 0. 1.
               1. 0. 1.
               0. 1. 1.];
  case 'PR15',
    Xcorel1 = [0.  0.  -1.
               1.  0.  -1.
               0.  1.  -1.
               0.  0.   1.
               1.  0.   1.
               0.  1.   1.
               0.5 0.  -1.
               0.5 0.5 -1.
               0.  0.5 -1.
               0.5 0.   1.
               0.5 0.5  1.
               0.  0.5  1.
               0.  0.   0.
               1.  0.   0.
               0.  1.   0.];
  case 'TET4',
    Xcorel1 = [0. 0. 0.
               1. 0. 0.
	       0. 1. 0.
	       0. 0. 1.];
  case 'TE10',
%   cf gmsh
    Xcorel1 = [0.  0.  0.
               1.  0.  0.
	       0.  1.  0.
	       0.  0.  1.
               0.5 0.  0.
               0.5 0.  0.5
               0.  0.  0.5
               0.  0.5 0.
               0.  0.5 0.5
               0.5 0.5 0.];
  case 'CUB8',
    Xcorel1 = [-1. -1. -1.
                1. -1. -1.
                1.  1. -1.
               -1.  1. -1.
               -1. -1.  1.
                1. -1.  1.
                1.  1.  1.
               -1.  1.  1.];
  case 'CU20',
    Xcorel1 = [-1. -1. -1.
                1. -1. -1.
                1.  1. -1.
               -1.  1. -1.
               -1. -1.  1.
                1. -1.  1.
                1.  1.  1.
               -1.  1.  1.
                0. -1. -1.
                1.  0. -1.
                0.  1. -1.
               -1.  0. -1.
                0. -1. 1.
                1.  0. 1.
                0.  1. 1.
               -1.  0. 1.
               -1. -1.  0.
                1. -1.  0.
                1.  1.  0.
               -1.  1.  0.];
  otherwise,
    type1
    error('Type of element not implemented yet')
end
