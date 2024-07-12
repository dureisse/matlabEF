function [PHI,DPHI,COOR,WEIGHT,varargout] = QUA8_SegmentIntegr3( ...
                                              support1,mode1,varargin);
% Integration information for QUA8 geometric element
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 09 / 2004
%
% Element fini = element geometrique (QUA8) + formulation + mode d'analyse
%
% [PHI,DPHI,COOR,WEIGHT] = QUA8_SegmentIntegr3(support1,mode1);
% [PHI,DPHI,COOR,WEIGHT] = QUA8_SegmentIntegr3(support1,mode1,intge2);
% [PHI,DPHI,COOR,WEIGHT,DDPHI] = QUA8_SegmentIntegr3(support1,mode1);
% [PHI,DPHI,COOR,WEIGHT,DDPHI] = QUA8_SegmentIntegr3(support1,mode1,intge2);
%
% Inputs
%   support1    Support d'integration
%               (RIGIDITE,MASSE) en formulation ELASTIQUE,
%               ou (COMPRESSIBILITE,PERMEABILITE) en formulation FLUIDE,
%               ou (CAPACITE,CONDUCTIVITE) en formulation THERMIQUE
%   mode1       Mode d'analyse (COPL,DEPL,TRID)
% Optional inputs
%   intge2      Segment d'integration a surcharger
% Outputs
%   PHI(nbnni,nbptg)            Valeur des fcts de base aux pts d'int.
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
X774 = .774596669241483D0;
P555 = .555555555555555555D0;
P888 = .888888888888888888D0;

if (strcmp(mode1,liste_mode{2}) || strcmp(mode1,liste_mode{3}) || ...
    strcmp(mode1,TRID))
%
%   COPL Plane stress or DEPL Plane strain or TRID Tridimensional
%   """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    idimr = 2; % nombre de coordonnees
    nbnni = 8; % Nombre de fonctions de base differentes
    if (strcmp(support1,liste_intg{1}) || ...
        strcmp(support1,liste_intg{4}) || ...
        strcmp(support1,liste_intg{7}))
%       RIGIDITE ou PERMEABILITE ou CONDUCTIVITE
        nbptg = 9; % Nombre de points d'integration
        COOR   = [-X774 -X774
	           0.   -X774
		   X774 -X774
                  -X774 0.
		   0.   0.
		   X774 0.
                  -X774 X774
		   0.   X774
		   X774 X774]; % Leurs coordonnees
        WEIGHT = [P555*P555
                  P888*P555
                  P555*P555
		  P888*P555
		  P888*P888
		  P888*P555
		  P555*P555
		  P888*P555
		  P555*P555]'; % Leur poids d'integration
    elseif (strcmp(support1,liste_intg{2}) || ...
            strcmp(support1,liste_intg{3}) || ...
            strcmp(support1,liste_intg{6}))
%       MASSE ou COMPRESSIBILITE ou CAPACITE
        disp(['Warning: underintegration for QUA8 and ' support1])
        nbptg = 9; % Cast3M integre avec ca... sous integration ?
        COOR   = [-X774 -X774
	           0.   -X774
		   X774 -X774
                  -X774 0.
		   0.   0.
		   X774 0.
                  -X774 X774
		   0.   X774
		   X774 X774]; % Leurs coordonnees
        WEIGHT = [P555*P555
                  P888*P555
                  P555*P555
		  P888*P555
		  P888*P888
		  P888*P555
		  P555*P555
		  P888*P555
		  P555*P555]'; % Leur poids d'integration
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
      PHI(:,ptg)    = EF_Shape('QUA8',COOR(ptg,:))'; % attention au '
      DPHI(:,:,ptg) = EF_Dshape('QUA8',COOR(ptg,:));
    end

else
    mode1
    error('mode not implemented')
end
