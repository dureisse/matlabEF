function [PHI,DPHI,COOR,WEIGHT,varargout] = PRI6_SegmentIntegr3( ...
                                             support1,mode1,varargin);
% Integration information for PRI6 geometric element
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 12 / 2006
%
% Retourne le segment d'integration pour l'element fini :
% element geometrique PRI6
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
X138 = 0.1381966011250105D0;
X585 = 0.5854101966249684D0;
X577 = 0.577350269189626D0;

if strcmp(mode1,'TRID')
%   TRID
    idimr = 3; % nombre de coordonnees
    nbnni = 6; % Nombre de fonctions de base
    if (strcmp(support1,liste_intg{1}) || ...
        strcmp(support1,liste_intg{4}) || ...
        strcmp(support1,liste_intg{7}))
%       RIGIDITE ou PERMEABILITE ou CONDUCTIVITE
        nbptg  = 6; % Nombre de points d'integration
        COOR   = [1./6. 1./6. -X577
                  2./3. 1./6. -X577
                  1./6. 2./3. -X577
                  1./6. 1./6. X577
                  2./3. 1./6. X577
                  1./6. 2./3. X577];
        WEIGHT = [1./6. 1./6. 1./6. 1./6. 1./6. 1./6.];
    elseif (strcmp(support1,liste_intg{2}) || ...
            strcmp(support1,liste_intg{3}) || ...
            strcmp(support1,liste_intg{6}))
%       MASSE ou COMPRESSIBILITE ou CAPACITE
        nbptg  = 6; % Nombre de points d'integration
        COOR   = [1./6. 1./6. -X577
                  2./3. 1./6. -X577
                  1./6. 2./3. -X577
                  1./6. 1./6. X577
                  2./3. 1./6. X577
                  1./6. 2./3. X577];
        WEIGHT = [1./6. 1./6. 1./6. 1./6. 1./6. 1./6.];
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
      PHI(:,ptg)    = EF_Shape('PRI6',COOR(ptg,:))'; % attention au '
      DPHI(:,:,ptg) = EF_Dshape('PRI6',COOR(ptg,:));
    end

else
  mode1
  error('mode not implemented')
end
