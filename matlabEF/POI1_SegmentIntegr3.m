function [PHI,DPHI,COOR,WEIGHT] = POI1_SegmentIntegr3(mode1);
% Integration information for POI1 element
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 21 / 01 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 05 / 02 / 2003
%   On ote la notion de dimension dans l'espace de reference (idimr)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 26 / 04 / 2004
%   Passage au 3D TRID
%
% Retourne le segment d'integration pour l'element fini :
% element geometrique POI1
% Independant de la formulation, du support d'integration :
% on somme de facon discrete au point considere.
%
% Entree
% """"""
% mode1 : mode d'analyse (COPL,DEPL)

GlobalVar;

if strcmp(mode1,liste_mode{2}) | strcmp(mode1,liste_mode{3}) | ...
   strcmp(mode1,liste_mode{4}) | strcmp(mode1,liste_mode{1})
%   COPL Plane stress or DEPL Plane strain or BARR Bar element
%   or TRID
%%    idimr = 2; % nombre de coordonnees
    idimr = 0; % nombre de coordonnees
    nbnni = 0; % Nombre de fonctions de base
    nbptg = 1; % Nombre de points d'integration
    COOR   = [0. 0.];
    WEIGHT = [1];
    PHI    = zeros(nbnni,nbptg);
    DPHI   = zeros(idimr,nbnni,nbptg);
else
    mode1
    error('mode not implemented')
end
