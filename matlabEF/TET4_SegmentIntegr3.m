function [PHI,DPHI,COOR,WEIGHT,varargout] = TET4_SegmentIntegr3( ...
                                             support1,mode1,varargin);
% Integration information for TET4 geometric element
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 24 / 08 / 2006
%
% Retourne le segment d'integration pour l'element fini :
% element geometrique TET4
% formulation ELASTIQUE ou FLUIDE ou THERMIQUE
% support d'integration support1 (RIGIDITE,MASSE) ou
% (COMPRESSIBILITE,PERMEABILITE) ou
% (CAPACITE,CONDUCTIVITE)
% mode d'analyse mode1 (TRID)

if nargin==2
  clear intge2;
elseif nargin==3
  intge2 = varargin{1};
else
  error('Bad number of arguments')
end
nout = nargout-4;

GlobalVar;

if strcmp(mode1,'TRID')
%   TRID
    idimr = 3; % nombre de coordonnees
    nbnni = 4; % Nombre de fonctions de base
    if (strcmp(support1,liste_intg{1}) || ...
        strcmp(support1,liste_intg{4}) || ...
        strcmp(support1,liste_intg{7}))
%       RIGIDITE ou PERMEABILITE ou CONDUCTIVITE
        nbptg  = 1; % Nombre de points d'integration
        [COOR,WEIGHT] = TET_Integr(nbptg);
    elseif (strcmp(support1,liste_intg{2}) || ...
            strcmp(support1,liste_intg{3}) || ...
            strcmp(support1,liste_intg{6}))
%       MASSE ou COMPRESSIBILITE ou CAPACITE
        nbptg  = 4;
        [COOR,WEIGHT] = TET_Integr(nbptg);
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
      PHI(:,ptg)    = EF_Shape('TET4',COOR(ptg,:))'; % attention au '
      DPHI(:,:,ptg) = EF_Dshape('TET4',COOR(ptg,:));
    end

else
  mode1
  error('mode not implemented')
end
