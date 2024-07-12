function [PHI,DPHI,COOR,WEIGHT,varargout] = TRI3_SegmentIntegr3( ...
                                             support1,mode1,varargin);
% Integration information for TRI3 geometric element
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 30 / 07 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 11 / 2002
%   Evite duplication des fcts de base
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 23 / 01 / 2003
%   Utilisation de EF_Shape, EF_Dshape
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 02 / 04 / 2003
%   Ajout du support COMPRESSIBILITE
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 05 / 2003
%   Possibilite de surcharger les points et poids d'integration
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 23 / 09 / 2003
%   Ajout des supports CAPACITE et CONDUCTIVITE
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 07 / 2004
%   Ajout du mode d'analyse DKIR
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 31 / 07 / 2007
%   Ajout du support ADVECTION
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 09 / 2007
%   Ajout du mode AXIS
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 11 / 2007
%   Factorisation de la definition des points d'integration
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 05 / 2008
%   Ajout du support COMPRESSIBILITE en DKT (pour ChamToChmno)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 07 / 2008
%   Ajout du mode TRID
% DUREISSEIX David  LaMCoS                           le 08 / 10 / 2010
%   Ajout du mode d'analyse DSHE
%
% Retourne le segment d'integration pour l'element fini :
% element geometrique TRI3
% formulation ELASTIQUE ou FLUIDE ou THERMIQUE ou TRANSPORT
% support d'integration support1 (RIGIDITE,MASSE) ou
% (COMPRESSIBILITE,PERMEABILITE) ou
% (CAPACITE,CONDUCTIVITE) ou
% MASSE,ADVECTION,FLUX
% mode d'analyse mode1 (COPL,DEPL,DKIR,DSHE,AXIS)
% (tout ne marche pas pour tout...)
%
% [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1);
% [PHI,DPHI,COOR,WEIGHT,DDPHI] = TRI3_SegmentIntegr3(support1,mode1);
% [PHI,DPHI,COOR,WEIGHT,***] = TRI3_SegmentIntegr3('FLUX',mode1);

if nargin==2
  clear intge2;
elseif nargin==3
  intge2 = varargin{1};
else
  error('Bad number of arguments')
end
nout = nargout-4;

GlobalVar;
%
%%%     INTEGRATION 7 POINTS TRIANGLE
%%      X101=.1012865073235D0;
%%      X797=.7974269853531D0;
%%      X470=.4701420641051D0;
%%      X059=.0597158717898D0;
%%      P125=.1259391805448D0;
%%      P132=.1323941527885D0;

if (strcmp(mode1,'COPL') || strcmp(mode1,'DEPL') || ...
    strcmp(mode1,'TRID'))
%
% COPL Plane stress, DEPL Plane strain or TRID Tridimensional
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  idimr = 2; % nombre de coordonnees
  nbnni = 3; % Nombre de fonctions de base
  FELEM = 'TRI3'; % Type d'element fini
  if (strcmp(support1,liste_intg{1}) || ...
      strcmp(support1,liste_intg{4}) || ...
      strcmp(support1,liste_intg{7}))
%   RIGIDITE ou PERMEABILITE ou CONDUCTIVITE
    nbptg  = 1; % Nombre de points d'integration
    [COOR,WEIGHT] = TRI_Integr(nbptg);
  elseif (strcmp(support1,liste_intg{2}) || ...
          strcmp(support1,liste_intg{3}) || ...
          strcmp(support1,liste_intg{6}))
%   MASSE ou COMPRESSIBILITE ou CAPACITE
    nbptg  = 3; % cf Zienkiewicz
    [COOR,WEIGHT] = TRI_Integr(nbptg);
  elseif strcmp(support1,liste_intg{8})
%   ADVECTION
    nbptg  = 1; % Like for Finite Volumes: centroid
    [COOR,WEIGHT] = TRI_Integr(nbptg);
  elseif strcmp(support1,liste_intg{9})
%   FLUX (on the boundary, edges numbered as for Gmsh)
    disp('Warning: TRANSPORT as formulation, FLUX as support:')
    disp('  integration along the edges of the elements')
    nbptg  = 3;
    COOR   = [0.5 0.
              0.5 0.5
              0.  0.5];
    WEIGHT = [1. sqrt(2) 1.]; % Beware: weight for the edges!
  else
    support1
    error('support not implemented')
  end

elseif strcmp(mode1,'AXIS')

% AXIS Axisymmetric
% """""""""""""""""
  idimr = 2; % nombre de coordonnees
  nbnni = 3; % Nombre de fonctions de base
  FELEM = 'TRI3'; % Type d'element fini
  if (strcmp(support1,liste_intg{1}) || ...
      strcmp(support1,liste_intg{4}) || ...
      strcmp(support1,liste_intg{7}))
%   RIGIDITE ou PERMEABILITE ou CONDUCTIVITE
%    nbptg  = 4; % Nombre de points d'integration, cf Castem
%%   Version poids central negatif Zienckiewick p. 201
%%   A ne pas utiliser en non-lineaire
%    disp('TRI3_SegmentIntegr3 : Modifier integration negative')
%    COOR   = [1./3. 1./3.
%              1./5. 1./5.
%              3./5. 1./5.
%              1./5. 3./5.];
%    WEIGHT = [-27./96.  25./96.  25./96. 25./96.];
    disp('DD:TRI3_SegmentIntegr3: for AXIS, underintegration of RIGIDITY')
    nbptg  = 1; % Nombre de points d'integration, cf Castem
    [COOR,WEIGHT] = TRI_Integr(nbptg);

  elseif (strcmp(support1,liste_intg{2}) || ...
          strcmp(support1,liste_intg{3}) || ...
          strcmp(support1,liste_intg{6}))
%   MASSE ou COMPRESSIBILITE ou CAPACITE
    nbptg  = 7; % cf Castem et Bathe p. 280
    [COOR,WEIGHT] = TRI_Integr(nbptg);
  else
    support1
    error('support not implemented for AXIS')
  end

elseif strcmp(mode1,liste_mode{7})
%
% DKIR: Plane Discrete Kirchhoff plate (TRI3 for membrane, DKT for bending)
% """"""""""""""""""""""""""""""""""""
  idimr = 2; % Number of coordinates in reference (natural) space
  nbnni = 3+3; % Number of distinct shape functions
  FELEM = 'DKT'; % Type d'element fini
  if strcmp(support1,liste_intg{1})
%   RIGIDITE
    nbptg = 3; % Number of integration points (3 Hammer points, see Batoz)
    COOR   = [0.5 0.5
              0.  0.5
              0.5 0.]; % Coordinates of integration points (in ref. space)
%   cf Castem
    COOR   = [0.5 0.
              0.  0.5
              0.5 0.5]; % Coordinates of integration points (in ref. space)
    WEIGHT = [1./6. 1./6. 1./6.]; % Weights of integration points
  elseif strcmp(support1,liste_intg{5})
%   NOEUDS
    nbptg = 3; % Number of integration points (nodes)
    COOR   = [0. 0.
              1. 0.
              0. 1.]; % Coordinates of nodes (in ref. space)
    WEIGHT = [1./6. 1./6. 1./6.]; % Weights
  elseif strcmp(support1,'COMPRESSIBILITE')
    nbptg  = 3;
    [COOR,WEIGHT] = TRI_Integr(nbptg);
  else
    support1
    error('support not implemented bis')
  end

elseif strcmp(mode1,liste_mode{9})
%
% DSHE: Plane Discrete shear plate (TRI3 for membrane, DST for bending)
% """"""""""""""""""""""""""""""""
  idimr = 2; % Number of coordinates in reference (natural) space
  nbnni = 3+3; % Number of distinct shape functions
%                (3 corners as real notes, 3 edge centers as virtual nodes)
  FELEM = 'DST'; % Type d'element fini
  if strcmp(support1,'RIGIDITE')
    nbptg = 3; % Number of integration points
%   (3 Hammer points, see Cast3m, not the same order as in Batoz)
    COOR   = [0.5 0.
              0.  0.5
              0.5 0.5]; % Coordinates of integration points (in ref. space)
    WEIGHT = [1./6. 1./6. 1./6.]; % Weights of integration points
  elseif strcmp(support1,'NOEUDS')
    nbptg = 3; % Number of integration points (nodes)
    COOR   = [0. 0.
              1. 0.
              0. 1.]; % Coordinates of nodes (in ref. space)
    WEIGHT = [1./6. 1./6. 1./6.]; % Weights
  elseif strcmp(support1,'COMPRESSIBILITE')
    error('Not yet done...')
    nbptg  = 3;
    [COOR,WEIGHT] = TRI_Integr(nbptg);
  else
    support1
    error('support not implemented bis')
  end

else
  mode1
  error('mode not implemented')
end

% Surcharge eventuelle
  if exist('intge2')
    nbptg  = length(intge2.WEIGHT);
    COOR   = intge2.COOR;
    WEIGHT = intge2.WEIGHT;
  end

% Initialize
  PHI    = zeros(nbnni,nbptg);
  DPHI   = zeros(idimr,nbnni,nbptg);
  switch FELEM
    case {'DKT','DST'}
      DDPHI  = zeros(idimr,idimr,nbnni,nbptg);
    case 'TRI3',
    otherwise
      FELEM
      error('Finite element nor recognized')
  end
% Setup
  for ptg = 1:nbptg
    PHI(:,ptg) = EF_Shape(FELEM,COOR(ptg,:))'; % attention au '
    DPHI(:,:,ptg) = EF_Dshape(FELEM,COOR(ptg,:));
    if (strcmp(FELEM,'DKT') || strcmp(FELEM,'DST'))
      DDPHI(:,:,:,ptg) = EF_DDshape(FELEM,COOR(ptg,:));
    end
  end
  if nout == 1
    varargout(1) = {DDPHI};
  end
