function [PHI,DPHI,COOR,WEIGHT,varargout] = RAC2_SegmentIntegr3( ...
                                              support1,mode1,varargin);
% Integration information for RAC2 geometric element
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 05 / 2005
%
% Element fini = element geometrique (QUA4) + formulation + mode d'analyse
%
% [PHI,DPHI,COOR,WEIGHT] = RAC2_SegmentIntegr3(support1,mode1);
% [PHI,DPHI,COOR,WEIGHT] = RAC2_SegmentIntegr3(support1,mode1,intge2);
% [PHI,DPHI,COOR,WEIGHT,DDPHI] = RAC2_SegmentIntegr3(support1,mode1);
% [PHI,DPHI,COOR,WEIGHT,DDPHI] = RAC2_SegmentIntegr3(support1,mode1,intge2);
%
% Inputs
%   support1    Support d'integration (RIGIDITE) en formulation ELASTIQUE,
%   mode1       Mode d'analyse (JOIN)
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
O577=0.577350269189626;

if strcmp(mode1,liste_mode{8})
%
% JOIN Element joint 2D (epaisseur nulle)
% """""""""""""""""""""""""""""""""""""""
  idimr = 1; % Number of coordinates in reference (natural) space
  nbnni = 2; % Number of distinct shape functions
  if strcmp(support1,liste_intg{1})
%   RIGIDITE
    nbptg = 2;      % Number of integration points
    COOR   = 0.5*[1.-O577
                  1.+O577]; % Leurs coordonnees
    WEIGHT = [0.5 0.5]; % Leur poids d'integration
%    COOR   = [0.5]; % Coord. of integration pts (in ref. space)
%    WEIGHT = [1.];  % Weights of integration points
  elseif strcmp(support1,liste_intg{5}) 
%   NOEUDS
    nbptg = 2; % Number of integration points at nodes
    COOR   = [0.
              1.]; % Coordinates of integr. points (in ref. space)
    WEIGHT = 0.5*[1. 1.]; % Weights of integration points
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
%  DDPHI  = zeros(idimr,idimr,nbnni,nbptg);
  for ptg = 1:nbptg
    PHI(:,ptg)       = EF_Shape('RAC2',COOR(ptg,:))'; % attention au '
    DPHI(:,:,ptg)    = EF_Dshape('RAC2',COOR(ptg,:));
%    DDPHI(:,:,:,ptg) = EF_DDshape('RAC2',COOR(ptg,:));
  end
%  if nout == 1
%    varargout(1) = {DDPHI};
%  end
  
else
    mode1
    error('mode not implemented')
end
