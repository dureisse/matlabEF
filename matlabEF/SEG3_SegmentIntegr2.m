function [PHI,DPHI,COOR,WEIGHT] = SEG3_SegmentIntegr2(support1,mode1);
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 05 / 12 / 2002

% Retourne le segment d'integration pour l'element fini :
% element geometrique SEG3 plonge dans un espace plus grand
% (idim = 2 ou 3)
% formulation ELASTIQUE
% support d'integration support1 (RIGIDITE,MASSE)
% mode d'analyse mode1 (COPL,DEPL)

% Rappel : integration de Gauss en 1D
% p points de Gauss permettent d'integrer exactement des poly de degre
% 2p-1

GlobalVar;
O577=0.577350269189626;
X774=0.774596669241483;
P555=5./9.;
P888=8./9.;

if strcmp(mode1,liste_mode{2}) | strcmp(mode1,liste_mode{3})
%   COPL Plane stress or DEPL Plane strain
    idimr = 1; % nombre de coordonnees
    nbnni = 3; % Nombre de fonctions de base differentes
    if strcmp(support1,liste_intg{1})
%       RIGIDITE
        nbptg = 2; % Nombre de points d'integration
        COOR   = [-O577
                   O577]; % Leurs coordonnees
        WEIGHT = [1. 1.]; % Leur poids d'integration
    elseif strcmp(support1,liste_intg{2})
%       MASSE
        nbptg = 3;
        COOR   = [-X774
                   0.
                   X774]; % Leurs coordonnees
        WEIGHT = [P555 P888 P555]; % Leur poids d'integration
    else
	 support1
	 error('support not implemented')
    end
    PHI    = zeros(nbnni,nbptg);
    DPHI   = zeros(idimr,nbnni,nbptg);
    for ptg = 1:nbptg
      x = COOR(ptg,1);
      PHI(:,ptg) = [0.5*x*(x-1.)
                    0.5*x*(x+1)
                    (1.-x)*(1.+x)];
      DPHI(:,:,ptg) = [(x-0.5) (x+0.5) (-2.*x)];
    end
else
    mode1
    error('mode not implemented')
end
