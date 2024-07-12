function [PHI,DPHI,COOR,WEIGHT,varargout] = SEG2_SegmentIntegr3( ...
                                              support1,mode1,varargin);
% Integration information for 2-node segment geometric element
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 05 / 12 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 13 / 12 / 2002
%   ajout FLUIDE
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 05 / 02 / 2003
%   ajout mode BARR
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 21 / 07 / 2003
%   Possibilite de surcharger les points et poids d'integration
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 12 / 2003
%   ajout mode TIMO poutre Timoshenko, courbure negligee
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 12 / 2003
%   ajout mode POUT poutre Euler-Bernoulli, courbure negligee
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 19 / 06 / 2004
%   ajout du mode TRID
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 07 / 2007
%   Ajout du mode AXIS
% DUREISSEIX David  LaMCoS                           le 13 / 06 / 2014
%   Ajout THERMIQUE,CONDUCTIVITE
%
% Retourne le segment d'integration pour l'element fini :
% element geometrique SEG2 plonge dans un espace plus grand
% (idim = 2 ou 3)
% formulation ELASTIQUE ou FLUIDE ou THERMIQUE
%
% Entrees
%   support1	support d'integration (RIGIDITE,MASSE,NOEUDS) ou
%               (COMPRESSIBILITE,PERMEABILITE) ou (CAPACITE,CONDUCTIVITE)
%   mode1       mode d'analyse (cf ci-dessous)
% Entrees optionnelles
%   intge2	segment d'integration elementaire a surcharger
% Sorties
%   PHI(nbnni,nbptg)		Values of the nbnni shape functions
%                               at the nbptg integr. points
%   DPHI(idimr,nbnni,nbptg)	Values of the derivates (w.r.t. idimr
%                               coordinates) of the nbnni shape
%                               functions at the nbptg integr. points
%   COOR(nbptg,idimr)		Reference coordinates of integr. points
%   WEIGHT(nbptg)		Their integration weight
% Sorties optionnelles
%   DDPHI(idimr,idimr,nbnni,nbptg)
%				Values of the second derivates
%                               (w.r.t. idimr coordinates) of the nbnni
%                               shape functions at the nbptg integr. pts
%                               This is provided for HER2 finite element
%
% mode d'analyse mode1        element fini associe 
% COPL,DEPL,TRID,BARR,TIMO    SEG2
% POUT                        HER2 (Hermite 2 noeuds)
%
% nbptg : nombre de points d'integration
% nbnni : nombre de fonctions de base differentes


% Pour information : formules d'integration a 3 et 4 points de Gauss
% d'apres cast3m
%
%nbptg = 3;
%P555=.555555555555555555D0;
%P888=.888888888888888888D0;
%X774=.774596669241483D0;
%COOR = [-1.*X774 ; 0. ; X774];
%WEIGHT = [P555 P888 P555];
%
%nbptg = 4;
%X861=.861136311594053D0;
%X339=.339981043584856D0;
%P347=.347854845137454D0;
%P652=.652145154862546D0;
%COOR = [-1.*X861 ; -1.*X339 ; X339 ; X861];
%WEIGHT = [P347 P652 P652 P347];


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
O577=0.577350269189626;

idimr = 1; % nombre de coordonnees

if (strcmp(mode1,'COPL') || strcmp(mode1,'DEPL') || ...
    strcmp(mode1,'TRID') || strcmp(mode1,'BARR'))
%
% COPL Plane stress or DEPL Plane strain or TRID or BARR Bar element
% """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  nbnni = 2;      % Nombre de fonctions de base differentes
  FELEM = 'SEG2'; % Type d'element fini
  if (strcmp(support1,liste_intg{1}) || ...
      strcmp(support1,liste_intg{4}) || ...
      strcmp(support1,'CONDUCTIVITE'))
%   RIGIDITE ou PERMEABILITE ou CONDUCTIVITE
    nbptg = 1; % Nombre de points d'integration
    COOR   = [0.]; % Leurs coordonnees
    WEIGHT = [2.]; % Leur poids d'integration
  elseif (strcmp(support1,liste_intg{2}) || ...
          strcmp(support1,liste_intg{3}) || ... 
          strcmp(support1,liste_intg{6}))
%   MASSE ou COMPRESSIBILITE ou CAPACITE
    nbptg = 2;
    COOR   = [-O577
               O577]; % Leurs coordonnees
    WEIGHT = [1. 1.]; % Leur poids d'integration
  elseif strcmp(support1,liste_intg{5})
%   NOEUDS
    nbptg = 2;        % Number of integration points
    COOR  = [-1.
              1.];    % Coordinates in reference space
    WEIGHT = [1. 1.]; % Weights
  else
    support1
    error('support not implemented')
  end

elseif strcmp(mode1,'AXIS')
%
% AXIS Axisymmetry
% """"""""""""""""
  nbnni = 2; % Nombre de fonctions de base differentes
  FELEM = 'SEG2';
%  if (strcmp(support1,liste_intg{1}) | strcmp(support1,liste_intg{4}))
%%   RIGIDITE ou PERMEABILITE
%    nbptg = 1; % Nombre de points d'integration
%    COOR   = [0.]; % Leurs coordonnees
%    WEIGHT = [2.]; % Leur poids d'integration
  if (strcmp(support1,liste_intg{2}) || ...
      strcmp(support1,liste_intg{3}) || ... 
      strcmp(support1,liste_intg{6}))
%   MASSE ou COMPRESSIBILITE ou CAPACITE
    disp('Warning: SEG2 + AXIS reduced integration for MASS')
    nbptg = 2;
    COOR   = [-O577
               O577]; % Leurs coordonnees
    WEIGHT = [1. 1.]; % Leur poids d'integration
  elseif strcmp(support1,liste_intg{5})
%   NOEUDS
    nbptg = 2;        % Number of integration points
    COOR  = [-1.
              1.];    % Coordinates in reference space
    WEIGHT = [1. 1.]; % Weights
  else
    support1
    error('support not implemented - bis')
  end

elseif strcmp(mode1,'TIMO')
%
% TIMO Timoshenko beam (straight)
% """"""""""""""""""""
  nbnni = 2; % Nombre de fonctions de base differentes
  FELEM = 'SEG2';
  if strcmp(support1,liste_intg{1})
%   RIGIDITE
%   Integration reduite
    disp('Warning: SEG2 + TIMO reduced integration')
    nbptg = 1; % Nombre de points d'integration
    COOR   = [0.]; % Leurs coordonnees
    WEIGHT = [2.]; % Leur poids d'integration
%%   Integration exacte (2 ptg car presence de la rotation)
%    disp('Warning: SEG2 + TIMO exact integration: possible shear locking')
%    nbptg = 2;
%    COOR   = [-O577
%               O577]; % Leurs coordonnees
%    WEIGHT = [1. 1.]; % Leur poids d'integration
  elseif strcmp(support1,liste_intg{5})
%   NOEUDS
    nbptg = 2;        % Number of integration points
    COOR  = [-1.
              1.];    % Coordinates in reference space
    WEIGHT = [1. 1.]; % Weights
  else
    mode1
    support1
    error('support not implemented - ter')
  end

elseif strcmp(mode1,'POUT')
%
% POUT Euler-Bernoulli beam (straight): Hermite element
% """""""""""""""""""""""""
  nbnni = 6; % Nombre de fonctions de base differentes
%             (Hermite + transformation)
  FELEM = 'HER2';
  if strcmp(support1,liste_intg{1})
%   RIGIDITE
    nbptg = 2; % Nombre de points d'integration
    COOR   = [-O577
               O577]; % Leurs coordonnees
    WEIGHT = [1. 1.]; % Leur poids d'integration
  elseif strcmp(support1,'MASSE')
%   MASSE
%   Integration reduite
    disp('Warning: SEG2 + POUT reduced integration for MASS')
    nbptg = 2; % Nombre de points d'integration
    COOR   = [-O577
               O577]; % Leurs coordonnees
    WEIGHT = [1. 1.]; % Leur poids d'integration
  elseif strcmp(support1,liste_intg{5})
%   NOEUDS
    nbptg = 2;        % Number of integration points
    COOR  = [-1.
              1.];    % Coordinates in reference space
    WEIGHT = [1. 1.]; % Weights
  else
    support1
    error('support not implemented for POUT')
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

% Values of the nbnni shape functions at the nbptg integr. points
PHI    = zeros(nbnni,nbptg);
% Values of the derivates (w.r.t. idimr coordinates) of the 
% nbnni shape functions at the nbptg integr. points
DPHI   = zeros(idimr,nbnni,nbptg);

switch FELEM
  case 'SEG2',
    for ptg = 1:nbptg
      x = COOR(ptg,1);
      PHI(:,ptg) = 0.5 * [(1.-x)
                          (1.+x)];
      DPHI(:,:,ptg) = 0.5 * [-1. 1.];
    end

  case 'HER2',
%   Values of the second derivates (w.r.t. idimr coordinates) of the 
%   nbnni shape functions at the nbptg integr. points
    DDPHI = zeros(idimr,idimr,nbnni,nbptg);
    for ptg = 1:nbptg
      x = COOR(ptg,1);
      PHI(:,ptg) = [0.25*(2.+x)*((x-1.)^2)
                    0.25*(2.-x)*((x+1.)^2)
                    0.25*(1.+x)*((x-1.)^2)
                    0.25*(x-1.)*((x+1.)^2)
                    0.5*(1.-x)
                    0.5*(1.+x)];
      DPHI(:,:,ptg) = [0.75*(x-1.)*(x+1.) ...
                       0.75*(x+1.)*(1.-x) ...
                       0.25*(x-1.)*(3.*x+1.)  ...
                       0.25*(x+1.)*(3.*x-1.) ...
                      -0.5 ...
                       0.5];
      DDPHI(1,1,:,ptg) = 1.5 * [x -x (x-1./3.) (x+1./3.) 0. 0.];
    end
    if nout == 1
      varargout(1) = {DDPHI};
    end

  otherwise,
    FELEM
    error('Unknown finite element')
end
