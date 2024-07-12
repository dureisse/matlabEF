function [PHI,DPHI,COOR,WEIGHT,varargout] = CUB8_SegmentIntegr3( ...
                                              support1,mode1,varargin);
% Integration information for CUB8 geometric element
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 13 / 12 / 2004
%
% Element fini = element geometrique (CUB8) + formulation + mode d'analyse
%
% [PHI,DPHI,COOR,WEIGHT] = CUB8_SegmentIntegr3(support1,mode1);
% [PHI,DPHI,COOR,WEIGHT] = CUB8_SegmentIntegr3(support1,mode1,intge2);
% [PHI,DPHI,COOR,WEIGHT,DDPHI] = CUB8_SegmentIntegr3(support1,mode1);
% [PHI,DPHI,COOR,WEIGHT,DDPHI] = CUB8_SegmentIntegr3(support1,mode1,intge2);
%
% Inputs
%   support1    Support d'integration (RIGIDITE,MASSE) en formulation ELASTIQUE,
%               ou (COMPRESSIBILITE,PERMEABILITE) en formulation FLUIDE,
%               ou (CAPACITE,CONDUCTIVITE) en formulation THERMIQUE
%   mode1       Mode d'analyse (TRID)
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
  nargin
  error('Bad number of arguments')
end
nout = nargout-4;

GlobalVar;
X577 = 0.577350269189626;

if strcmp(mode1,liste_mode{1})
%   TRID Tridimensionnel
    idimr = 3; % nombre de coordonnees
    nbnni = 8; % Nombre de fonctions de base differentes
    if (strcmp(support1,liste_intg{1}) || ...
        strcmp(support1,liste_intg{4}) || ...
        strcmp(support1,liste_intg{7}))
%       RIGIDITE ou PERMEABILITE ou CONDUCTIVITE
        nbptg = 8; % Nombre de points d'integration
        COOR   = [-X577 -X577 -X577
                   X577 -X577 -X577
                   X577  X577 -X577
                  -X577  X577 -X577
                  -X577 -X577  X577
                   X577 -X577  X577
                   X577  X577  X577
                  -X577  X577  X577]; % Leurs coordonnees
        WEIGHT = [1. 1. 1. 1. 1. 1. 1. 1.]; % Leur poids d'integration
    elseif (strcmp(support1,liste_intg{2}) || ...
            strcmp(support1,liste_intg{3}) || ...
            strcmp(support1,liste_intg{6}))
%       MASSE ou COMPRESSIBILITE ou CAPACITE
        disp(['Warning: underintegration for CUB8 and ' support1])
        nbptg = 8; % Cast3M integre avec ca... sous integration ?
        COOR   = [-X577 -X577 -X577
                   X577 -X577 -X577
                   X577  X577 -X577
                  -X577  X577 -X577
                  -X577 -X577  X577
                   X577 -X577  X577
                   X577  X577  X577
                  -X577  X577  X577]; % Leurs coordonnees
        WEIGHT = [1. 1. 1. 1. 1. 1. 1. 1.]; % Leur poids d'integration
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
      x = COOR(ptg,1); y = COOR(ptg,2);
      PHI(:,ptg) = EF_Shape('CUB8',COOR(ptg,:))'; % attention au '
      DPHI(:,:,ptg) = EF_Dshape('CUB8',COOR(ptg,:));
    end

else
    mode1
    error('mode not implemented')
end
