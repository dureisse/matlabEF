function [PHI,DPHI,COOR,WEIGHT] = CU20_SegmentIntegr3( ...
                                             support1,mode1,varargin);
% Integration information for CU20 geometric element
%
% DUREISSEIX David  LaMCoS                         le 02 / 10 / 2019
%
% Inputs
%   support1    Integration support
%   mode1       Analysis mode
% Optional input
%   intge2      Integration segment to be overloaded
% Outputs
%    PHI(nbnni,nbptg)        Values of shape functions at integration nodes
%    DPHI(idimr,nbnni,nbptg) Values of shape function derivatives
%    COOR(nbptg,idimr)       Reference coordinates of integration points
%    WEIGHT(nbptg)           Weights of integration points
%
% Retourne le segment d'integration pour l'element fini :
% element geometrique CU20
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

% cf Cast3M donred.eso
X774 = 0.774596669241483;
P555 = 0.555555555555555555;
P888 = 0.888888888888888888;
 

if strcmp(mode1,'TRID')
%   TRID
    idimr = 3; % nombre de coordonnees
    nbnni = 20; % Nombre de fonctions de base
    if (strcmp(support1,'RIGIDITE') || strcmp(support1,'PERMEABILITE') || ...
        strcmp(support1,'CONDUCTIVITE'))
%     RIGIDITE ou PERMEABILITE ou CONDUCTIVITE
      nbptg  = 27; % cf Cast3m elquoi.eso, donred.eso
      % Integration point coordinates
      COOR = zeros(nbptg,idimr); % QSIGAU = COOR(:,1); ETAGAU = COOR(:,2); DZEGAU = COOR(:,3);
      for IA = 1:3:25
          COOR(IA,1)   = -X774;
          COOR(IA+2,1) = X774;
      end
      for IA = 1:9
          COOR(IA,3)    = -X774;
          COOR(IA+18,3) = X774;
      end
      for IA = 1:3
          COOR(IA,2)    = -X774;
          COOR(IA+6,2)  = X774;
          COOR(IA+9,2)  = -X774;
          COOR(IA+15,2) = X774;
          COOR(IA+18,2) = -X774;
          COOR(IA+24,2) = X774;
      end
      % Integration point weights
      WEIGHT(1) = P555*P555;
      WEIGHT(3) = WEIGHT(1);
      WEIGHT(7) = WEIGHT(1);
      WEIGHT(9) = WEIGHT(1);
      WEIGHT(2) = P888*P555;
      WEIGHT(4) = WEIGHT(2);
      WEIGHT(6) = WEIGHT(2);
      WEIGHT(8) = WEIGHT(2);
      WEIGHT(5) = P888*P888;
      for IA = 1:9
          XX = WEIGHT(IA);
          WEIGHT(IA)    = XX*P555;
          WEIGHT(IA+9)  = XX*P888;
          WEIGHT(IA+18) = WEIGHT(IA);
      end
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

    % cf cast3m shape3.eso
    PHI  = zeros(nbnni,nbptg);
    DPHI = zeros(idimr,nbnni,nbptg);
    for ptg = 1:nbptg
%      x = COOR(ptg,1); y = COOR(ptg,2); z = COOR(ptg,3);
      PHI(:,ptg)    = EF_Shape('CU20',COOR(ptg,:))'; % attention au '
      DPHI(:,:,ptg) = EF_Dshape('CU20',COOR(ptg,:));
    end

else
  mode1
  error('mode not implemented')
end
