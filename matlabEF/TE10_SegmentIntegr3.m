function [PHI,DPHI,COOR,WEIGHT] = TE10_SegmentIntegr3( ...
                                             support1,mode1,varargin);
% Integration information for TE10 geometric element
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 24 / 11 / 2007
%
% Inputs
%   support1		Integration support
%   mode1		Analysis mode
% Optional input
%   intge2		Integration segment to be overloaded
% Outputs
%    PHI(nbnni,nbptg)		Values of shape functions at integration nodes
%    DPHI(idimr,nbnni,nbptg)	Values of shape function derivatives
%    COOR(nbptg,idimr)		Reference coordinates of integration points
%    WEIGHT(nbptg)		Weights of integration points
%
% Retourne le segment d'integration pour l'element fini :
% element geometrique TE10
% formulation ELASTIQUE ou FLUIDE ou THERMIQUE
% support d'integration support1 (RIGIDITE)
% mode d'analyse mode1 (TRID)

if nargin==2
  clear intge2;
elseif nargin==3
  intge2 = varargin{1};
else
  error('Bad number of arguments')
end

GlobalVar;

if strcmp(mode1,'TRID')
%   TRID
    idimr = 3; % nombre de coordonnees
    nbnni = 10; % Nombre de fonctions de base
    if strcmp(support1,'RIGIDITE') | strcmp(support1,'PERMEABILITE') | ...
       strcmp(support1,'CONDUCTIVITE') 
%     RIGIDITE ou PERMEABILITE ou CONDUCTIVITE
      nbptg  = 4; % cf Castem
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
      PHI(1:10,ptg)     = EF_Shape('TE10',COOR(ptg,:))'; % attention au '
      PHI(11:14,ptg)    = EF_Shape('TET4',COOR(ptg,:))'; % attention au '
      DPHI(:,1:10,ptg)  = EF_Dshape('TE10',COOR(ptg,:));
      DPHI(:,11:14,ptg) = EF_Dshape('TET4',COOR(ptg,:));
    end

else
  mode1
  error('mode not implemented')
end