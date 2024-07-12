function [PHI,DPHI,COOR,WEIGHT,varargout] = QUA4_SegmentIntegr3( ...
                                              support1,mode1,varargin);
% Integration information for QUA4 geometric element
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 31 / 07 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 11 / 2002
%   Evite duplication des fcts de base
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 23 / 01 / 2003
%   Utilisation de EF_Shape, EF_Dshape
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 04 / 2004
%   Possibilite de surcharger les points et poids d'integration
%   Ajout des supports CAPACITE et CONDUCTIVITE
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 03 / 08 / 2004
%   Ajout du mode d'analyse DKIR
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 12 / 2004
%   Ajout du mode d'analyse TRID
%
% Element fini = element geometrique (QUA4) + formulation + mode d'analyse
%
% [PHI,DPHI,COOR,WEIGHT] = QUA4_SegmentIntegr3(support1,mode1);
% [PHI,DPHI,COOR,WEIGHT] = QUA4_SegmentIntegr3(support1,mode1,intge2);
% [PHI,DPHI,COOR,WEIGHT,DDPHI] = QUA4_SegmentIntegr3(support1,mode1);
% [PHI,DPHI,COOR,WEIGHT,DDPHI] = QUA4_SegmentIntegr3(support1,mode1,intge2);
%
% Inputs
%   support1    Support d'integration (RIGIDITE,MASSE) en formulation ELASTIQUE,
%               ou (COMPRESSIBILITE,PERMEABILITE) en formulation FLUIDE,
%               ou (CAPACITE,CONDUCTIVITE) en formulation THERMIQUE
%   mode1       Mode d'analyse (COPL,DEPL,DKIR,TRID)
% Optional inputs
%   intge2      Segment d'integration a surcharger
% Outputs
%   PHI(nbnni,nbptg)            Valeur des fcts de base aux points d'integration
%   DPHI(idimr,nbnni,nbptg)     Valeur des derivees premieres
% Optional outputs
%   DDPHI(idimr,idimr,nbnni,nbptg)      Valeur des derivees secondes

if nargin==2
  clear intge2;
elseif nargin==3
  intge2 = varargin{1};
else
  error('Bad number of arguments')
end
nout = nargout-4;

GlobalVar;
r2 = sqrt(2.)/2.;
r3 = sqrt(3.)/3.;
r35 = sqrt(3./5.);

if strcmp(mode1,liste_mode{2}) || ...
   strcmp(mode1,liste_mode{3}) || ...
   strcmp(mode1,TRID)
%
%   COPL Plane stress or DEPL Plane strain or TRID Tridimensional
%   """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    idimr = 2; % nombre de coordonnees
    nbnni = 4; % Nombre de fonctions de base differentes
    if (strcmp(support1,liste_intg{1}) || ...
        strcmp(support1,liste_intg{4}) || ...
        strcmp(support1,liste_intg{7}))
%       RIGIDITE ou PERMEABILITE ou CONDUCTIVITE
        nbptg = 4; % Nombre de points d'integration
        COOR   = [-r3 -r3
                   r3 -r3
                   r3  r3
                  -r3  r3]; % Leurs coordonnees
        WEIGHT = [1. 1. 1. 1.]; % Leur poids d'integration
    elseif (strcmp(support1,liste_intg{2}) || ...
            strcmp(support1,liste_intg{3}) || ...
            strcmp(support1,liste_intg{6}))
%       MASSE ou COMPRESSIBILITE ou CAPACITE
        disp(['Warning: underintegration for QUA4 and ' support1])
        nbptg = 4; % Cast3M integre avec ca... sous integration ?
        COOR   = [-r3 -r3
                   r3 -r3
                   r3  r3
                  -r3  r3]; % Leurs coordonnees
        WEIGHT = [1. 1. 1. 1.]; % Leur poids d'integration
    else
	 support1
	 error('support not implemented')
    end

%   Surcharge eventuelle
    if exist('intge2')
      nbptg  = length(intge2.WEIGHT);
      COOR   = intge2.COOR;
      WEIGHT = intge2.WEIGHT;
    end

    PHI    = zeros(nbnni,nbptg);
    DPHI   = zeros(idimr,nbnni,nbptg);
    for ptg = 1:nbptg
%      x = COOR(ptg,1); y = COOR(ptg,2);
      PHI(:,ptg)    = EF_Shape('QUA4',COOR(ptg,:))'; % attention au '
      DPHI(:,:,ptg) = EF_Dshape('QUA4',COOR(ptg,:));
    end

elseif strcmp(mode1,liste_mode{7})
%
% DKIR Plane plate for Discrete Kirchhoff
% """""""""""""""""""""""""""""""""""""""
% (QUA4 for membrane, DKQ for bending)
  idimr = 2; % Number of coordinates in reference (natural) space
  nbnni = 4+4; % Number of distinct shape functions
  if strcmp(support1,liste_intg{1})
%   RIGIDITE
%    nbptg = 4; % Number of integration points (2x2 Gauss)
%    disp('Warning: underintegration for DKQ, rank should be OK')
%    COOR   = [-r3 -r3
%               r3 -r3
%               r3  r3
%              -r3  r3]; % Coordinates of integration points (in ref. space)
%    WEIGHT = [1. 1. 1. 1.]; % Weights of integration points
    nbptg = 9; % Number of integration points (3x3 Gauss)
    disp('Warning: high order integration for DKQ')
    COOR   = [-r35 -r35
               0.  -r35
               r35 -r35
              -r35  0.
               0.   0.
               r35  0.
              -r35  r35
               0.   r35
               r35  r35]; % Coordinates of integration points (in ref. space)
    WEIGHT = [25./81. 40./81. 25./81. ...
              40./81. 64./81  40./81. ...
              25./81. 40./81. 25./81.]; % Weights of integration points
  elseif strcmp(support1,liste_intg{5}) 
%   NOEUDS
    nbptg = 4; % Number of integration points (number of nodes)
    COOR   = [0. 0.
              1. 0.
              1. 1.
              0. 1.]; % Coordinates of integration points (in ref. space)
    WEIGHT = [1. 1. 1. 1.]; % Weights of integration points
  elseif strcmp(support1(1:7),'REGULAR')
%   REGULARn : sous-decoupage regulier avec nxn points d'integration
    n = str2num(support1(8:end));
    nbptg = n^2;
    disp(['Regular integration for DKQ, with ' int2str(nbptg) ' points'])
    COOR = zeros(nbptg,2);
    for i = 1:n
      xi = -1.+(1./n) + (i-1)*2./n;
      for j = 1:n
        yj = -1.+(1./n) + (j-1)*2./n;
        ptg = (i-1)*n+j;
        COOR(ptg,:) = [xi yj];
      end
    end
    WEIGHT = (4./nbptg) * ones(1,nbptg);
  elseif strcmp(support1(1:8),'RIGIDITE')
%   RIGIDITEn : sous-decoupage hierarchique avec nxn points d'integr.
    n = str2num(support1(9:end));
    nbptg = 4*(n^2);
    disp(['Hierarchical integration for DKQ rigidity, with ' ...
           int2str(nbptg) ' points'])
    coorr = [-r3 -r3
              r3 -r3
              r3  r3
             -r3  r3];
    COOR = zeros(nbptg,2);
    for i = 1:n
      xi = -1.+(1./n) + (i-1)*2./n;
      for j = 1:n
        yj = -1.+(1./n) + (j-1)*2./n;
        for iptg = 1:4
          ptg = ((i-1)*n+j-1)*4 + iptg;
          COOR(ptg,:) = ((1. / n) * coorr(iptg,:)) + [xi yj];
        end
      end
    end
    WEIGHT = (4./nbptg) * ones(1,nbptg);
  else
    support1
    error('support not implemented 2')
  end

% Surcharge eventuelle
  if exist('intge2')
    nbptg  = length(intge2.WEIGHT);
    COOR   = intge2.COOR;
    WEIGHT = intge2.WEIGHT;
  end

  PHI    = zeros(nbnni,nbptg);
  DPHI   = zeros(idimr,nbnni,nbptg);
  DDPHI  = zeros(idimr,idimr,nbnni,nbptg);
  for ptg = 1:nbptg
%    x = COOR(ptg,1); y = COOR(ptg,2);
    PHI(:,ptg)       = EF_Shape('DKQ',COOR(ptg,:))'; % attention au '
    DPHI(:,:,ptg)    = EF_Dshape('DKQ',COOR(ptg,:));
    DDPHI(:,:,:,ptg) = EF_DDshape('DKQ',COOR(ptg,:));
  end
  if nout == 1
    varargout(1) = {DDPHI};
  end
  
else
    mode1
    error('mode not implemented')
end
