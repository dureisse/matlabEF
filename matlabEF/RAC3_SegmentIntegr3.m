function [PHI,DPHI,COOR,WEIGHT,varargout] = RAC3_SegmentIntegr3( ...
                                              support1,mode1,varargin);
% Integration information for RAC3 geometric element
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 09 / 2005
%
% Element fini = element geometrique (QUA8) + formulation + mode d'analyse
%
% [PHI,DPHI,COOR,WEIGHT] = RAC3_SegmentIntegr3(support1,mode1);
% [PHI,DPHI,COOR,WEIGHT] = RAC3_SegmentIntegr3(support1,mode1,intge2);
% [PHI,DPHI,COOR,WEIGHT,DDPHI] = RAC3_SegmentIntegr3(support1,mode1);
% [PHI,DPHI,COOR,WEIGHT,DDPHI] = RAC3_SegmentIntegr3(support1,mode1,intge2);
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
X774=0.774596669241483;
P555=5./9.;
P888=8./9.;

if strcmp(mode1,liste_mode{8})
%
% JOIN Element joint 2D (epaisseur nulle)
% """""""""""""""""""""""""""""""""""""""
  idimr = 1; % Number of coordinates in reference (natural) space
  nbnni = 3; % Number of distinct shape functions
  if strcmp(support1,liste_intg{1})
%   RIGIDITE
    nbptg = 3;      % Number of integration points
    COOR   = [-X774
               0.
	       X774]; % Coordinates of integr. points (in ref. space)
    WEIGHT = [P555 P888 P555]; % Weights of integration points
  elseif strcmp(support1,liste_intg{5}) 
%   NOEUDS
    nbptg = 3; % Number of integration points at nodes = nb of nodes
    COOR   = [-1.
               0.
               1.]; % Coordinates of integr. points (in ref. space)
    WEIGHT = (2./3.)*[1. 1. 1.]; % Weights of integration points
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
    PHI(:,ptg)       = EF_Shape('RAC3',COOR(ptg,:))'; % attention au '
    DPHI(:,:,ptg)    = EF_Dshape('RAC3',COOR(ptg,:));
%    DDPHI(:,:,:,ptg) = EF_DDshape('RAC3',COOR(ptg,:));
  end
%  if nout == 1
%    varargout(1) = {DDPHI};
%  end
  
else
    mode1
    error('mode not implemented')
end
