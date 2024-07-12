function [PHI,DPHI,COOR,WEIGHT] = TET5b_Poreux_SegmentIntegr3( ...
                                             support1,mode1,varargin);
% Integration information for TET5b geometric element for porous case
%
% TET5b is a TET4 element plus a bubble-like 5th node at the centroid
% The "bubble" is a pyramid-like shape function (linear on each
% sub-triangle), non differentiable on the radial lines (!)
% Original shape functions of the TET4 element are modified to 
% preserve interpolation property, see Hughes (The FEM, linear
% static and dynamic FEA).
% Integration is selected as 1 integration point on each sub-triangle.
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
% element geometrique TET4
% formulation POREUX
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
    nbnni = 5+4; % Nombre de fonctions de base
    if strcmp(support1,'RIGIDITE')
%       RIGIDITE
        nbptg  = 4; % Nombre de points d'integration
%%        [COOR,WEIGHT] = TET_Integr(nbptg);
%       Centroid is 0.25 0.25 0.25
%       Au centre des sous-triangles
        COOR = [0.3125    0.3125    0.3125
                0.0625    0.3125    0.3125
                0.3125    0.0625    0.3125
                0.3125    0.3125    0.0625];
        WEIGHT = 1./24.*ones(1,4);
%       Integration point number i is the farthest from node i
%       Sub-triangle i contains integration point i
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
      PHI(1:5,ptg)    = EF_Shape('TET5b',COOR(ptg,:))'; % attention au '
      PHI(6:9,ptg)    = EF_Shape('TET4',COOR(ptg,:))'; % attention au '
      DPHI(:,1:5,ptg) = EF_Dshape('TET5b',COOR(ptg,:));
      DPHI(:,6:9,ptg) = EF_Dshape('TET4',COOR(ptg,:));
    end

else
  mode1
  error('mode not implemented')
end
