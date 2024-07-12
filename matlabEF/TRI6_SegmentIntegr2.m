function [PHI,DPHI,COOR,WEIGHT] = TRI6_SegmentIntegr2(support1,mode1,varargin);
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 11 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 25 / 11 / 2003
%   Ajout du support CAPACITE
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 27 / 05 / 2004
%   Possibilite de surcharger les points et poids d'integration
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 09 / 07 / 2007
%   Ajout mode AXIS
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 11 / 2007
%   Factorisation de la definition des points d'integration
%
% Retourne le segment d'integration pour l'element fini :
% element geometrique TRI6 isoparametrique
% (voir aussi TRI6d_Poreux_SegmentIntegr2)
% formulation ELASTIQUE
% support d'integration support1 (RIGIDITE) ou
% (CAPACITE)
% mode d'analyse mode1 (COPL,DEPL,AXIS)
%
%    3
%   | \
%   |  \
%   6   5
%   |    \
%   |     \
%   1__4___2
%

if nargin==2
  clear intge2;
elseif nargin==3
  intge2 = varargin{1};
else
  error('Bad number of arguments')
end

GlobalVar;

if (strcmp(mode1,liste_mode{2}) || strcmp(mode1,liste_mode{3}))
% COPL Plane stress or DEPL Plane strain
  idimr = 2; % nombre de coordonnees
  nbnni = 6; % Nombre de fonctions de base differentes

  if (strcmp(support1,liste_intg{1}) || ...
      strcmp(support1,liste_intg{4}) || ...
      strcmp(support1,liste_intg{7}))
%   RIGIDITE ou PERMEABILITE ou CONDUCTIVITE

if (1 == 0)
    nbptg = 4; % Nombre de points d'integration
%   Version poids central negatif Zienckiewick p. 201
%   A ne pas utiliser en non-lineaire
    disp('Modifier integration negative')
    COOR   = [1./3. 1./3.
              1./5. 1./5.
              3./5. 1./5.
              1./5. 3./5.];
    WEIGHT = [-27./96.  25./96.  25./96. 25./96.];
end

%   Version 3 points de Gauss a l'interieur de l'element
%   Il existe aussi une version 3 points d'integration au milieu
%   des aretes (ordre a confirmer).
    nbptg  = 3; % cf Zienkiewicz
    [COOR,WEIGHT] = TRI_Integr(nbptg);
%%    COOR   = [1./6. 1./6.
%%              2./3. 1./6.
%%              1./6. 2./3.];
%%    WEIGHT = [1./6. 1./6. 1./6.];

  elseif (strcmp(support1,liste_intg{2}) || ...
          strcmp(support1,liste_intg{3}) || ...
          strcmp(support1,liste_intg{6}))
%   MASSE ou COMPRESSIBILITE ou CAPACITE

    nbptg = 6; % Nombre de points d'integration
    [COOR,WEIGHT] = TRI_Integr(nbptg);
%%%   dans BATOZ tome 1 p200 : integration de poly d'ordre 4
%%    A = 0.445948490915965;
%%    B = 0.091576213509771;
%%    C = 0.111690794839005;
%%    D = 0.054975871827661;
%%    COOR = [ ...
%%      A          A
%%      (1. -2.*A) A
%%      A          (1.-2.*A)
%%      B          B
%%      (1.-2.*B)  B
%%      B          (1.-2.*B)];
%%    WEIGHT = [C C C D D D];

if (1 == 0)
    nbptg = 7; % Nombre de points d'integration
%   Version programmee dans Cast3M en thermique,
%   dans BATOZ tome 1 p200 : integration de poly d'ordre 5
    [COOR,WEIGHT] = TRI_Integr(nbptg);
%%    COOR = [...
%%      1./3.                 1./3.
%%      (9.-2.*sqrt(15.))/21. (6.+sqrt(15.))/21.
%%      (6.+sqrt(15.))/21.    (9.-2.*sqrt(15.))/21.
%%      (6.+sqrt(15.))/21.    (6.+sqrt(15.))/21.
%%      (9.+2.*sqrt(15.))/21. (6.-sqrt(15.))/21.
%%      (6.-sqrt(15.))/21.    (9.+2.*sqrt(15.))/21.
%%      (6.-sqrt(15.))/21.    (6.-sqrt(15.))/21.];
%%    WEIGHT = [ 9./80. ...
%%      (155.+sqrt(15.))/2400. (155.+sqrt(15.))/2400. (155.+sqrt(15.))/2400. ...
%%      (155.-sqrt(15.))/2400. (155.-sqrt(15.))/2400. (155.-sqrt(15.))/2400.];
end

  else
    support1
    error('support not implemented')
  end

% Surcharge eventuelle
  if exist('intge2')
    nbptg  = length(intge2.WEIGHT);
    COOR   = intge2.COOR;
    WEIGHT = intge2.WEIGHT;
  end

  PHI    = zeros(nbnni,nbptg);
  DPHI   = zeros(idimr,nbnni,nbptg);
  for ptg = 1:nbptg
    x = COOR(ptg,1); y = COOR(ptg,2);
    PHI(:,ptg) = [(1.-x-y)*(1.-2.*x-2.*y)
                  x*(2.*x-1.)
                  y*(2.*y-1.)
                  4.*x*(1.-x-y)
                  4.*x*y
                  4.*y*(1.-x-y)];
    DPHI(:,:,ptg) = [(4.*x+4.*y-3)  (4.*x+4.*y-3)
                     (4.*x-1.)      0.
                     0.             (4.*y-1.)
                     (4.-8.*x-4.*y) (-4.*x)
                     (4.*y)         (4.*x)
                     (-4.*y)        (4.-4.*x-8.*y)]';
  end

elseif strcmp(mode1,'AXIS')
% AXIS Axisymmetric
  idimr = 2; % nombre de coordonnees
  nbnni = 6; % Nombre de fonctions de base differentes

  if (strcmp(support1,liste_intg{1}) || ...
      strcmp(support1,liste_intg{4}) || ...
      strcmp(support1,liste_intg{7}))
%   RIGIDITE ou PERMEABILITE ou CONDUCTIVITE

    nbptg = 7; % Nombre de points d'integration
    [COOR,WEIGHT] = TRI_Integr(nbptg);
%%%   Version programmee dans Cast3M en thermique,
%%%   dans BATOZ tome 1 p200 : integration de poly d'ordre 5
%%    COOR = [...
%%      1./3.                 1./3.
%%      (9.-2.*sqrt(15.))/21. (6.+sqrt(15.))/21.
%%      (6.+sqrt(15.))/21.    (9.-2.*sqrt(15.))/21.
%%      (6.+sqrt(15.))/21.    (6.+sqrt(15.))/21.
%%      (9.+2.*sqrt(15.))/21. (6.-sqrt(15.))/21.
%%      (6.-sqrt(15.))/21.    (9.+2.*sqrt(15.))/21.
%%      (6.-sqrt(15.))/21.    (6.-sqrt(15.))/21.];
%%    WEIGHT = [ 9./80. ...
%%      (155.+sqrt(15.))/2400. (155.+sqrt(15.))/2400. (155.+sqrt(15.))/2400. ...
%%      (155.-sqrt(15.))/2400. (155.-sqrt(15.))/2400. (155.-sqrt(15.))/2400.];

  elseif (strcmp(support1,liste_intg{2}) || ...
          strcmp(support1,liste_intg{3}) || ...
          strcmp(support1,liste_intg{6}))
%   MASSE ou COMPRESSIBILITE ou CAPACITE

    nbptg = 7; % Nombre de points d'integration
    [COOR,WEIGHT] = TRI_Integr(nbptg);
%%%   Version programmee dans Cast3M en thermique,
%%%   dans BATOZ tome 1 p200 : integration de poly d'ordre 5
%%    COOR = [...
%%      1./3.                 1./3.
%%      (9.-2.*sqrt(15.))/21. (6.+sqrt(15.))/21.
%%      (6.+sqrt(15.))/21.    (9.-2.*sqrt(15.))/21.
%%      (6.+sqrt(15.))/21.    (6.+sqrt(15.))/21.
%%      (9.+2.*sqrt(15.))/21. (6.-sqrt(15.))/21.
%%      (6.-sqrt(15.))/21.    (9.+2.*sqrt(15.))/21.
%%      (6.-sqrt(15.))/21.    (6.-sqrt(15.))/21.];
%%    WEIGHT = [ 9./80. ...
%%      (155.+sqrt(15.))/2400. (155.+sqrt(15.))/2400. (155.+sqrt(15.))/2400. ...
%%      (155.-sqrt(15.))/2400. (155.-sqrt(15.))/2400. (155.-sqrt(15.))/2400.];

  else
    support1
    error('support not implemented')
  end

% Surcharge eventuelle
  if exist('intge2')
    nbptg  = length(intge2.WEIGHT);
    COOR   = intge2.COOR;
    WEIGHT = intge2.WEIGHT;
  end

  PHI    = zeros(nbnni,nbptg);
  DPHI   = zeros(idimr,nbnni,nbptg);
  for ptg = 1:nbptg
    x = COOR(ptg,1); y = COOR(ptg,2);
    PHI(:,ptg) = [(1.-x-y)*(1.-2.*x-2.*y)
                  x*(2.*x-1.)
                  y*(2.*y-1.)
                  4.*x*(1.-x-y)
                  4.*x*y
                  4.*y*(1.-x-y)];
    DPHI(:,:,ptg) = [(4.*x+4.*y-3)  (4.*x+4.*y-3)
                     (4.*x-1.)      0.
                     0.             (4.*y-1.)
                     (4.-8.*x-4.*y) (-4.*x)
                     (4.*y)         (4.*x)
                     (-4.*y)        (4.-4.*x-8.*y)]';
  end

else
  mode1
  error('mode not implemented')
end
