function [PHI,DPHI,COOR,WEIGHT] = TRI6d_Poreux_SegmentIntegr2(support1,mode1);
% DUREISSEIX David    LMGC EQUIPE SYSTEMES MULTICONTACTS  le 09 / 11 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 11 / 2002
%   Evite duplication des fcts de base
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 26 / 08 / 2003
%   Correction du nb de points de Gauss (4 points surabondants -> 3)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 09 / 2007
%   Ajout AXIS
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 11 / 2007
%   Factorisation de la definition des points d'integration

% Retourne le segment d'integration pour l'element fini :
% element geometrique TRI6 a bords droits
% (voir aussi TRI6_SegmentIntegr2)
% formulation POREUX
% support d'integration support1 (RIGIDITE)
% mode d'analyse mode1 (COPL,DEPL,AXIS)

GlobalVar;

%%%
%%%     INTEGRATION 7 POINTS TRIANGLE
%%      X101=.1012865073235D0;
%%      X797=.7974269853531D0;
%%      X470=.4701420641051D0;
%%      X059=.0597158717898D0;
%%      P125=.1259391805448D0;
%%      P132=.1323941527885D0;

if strcmp(mode1,'COPL') | strcmp(mode1,'DEPL')
% COPL Plane stress or DEPL Plane strain

  idimr  = 2;  % Nombre de coordonnees
  nbnnis = 6;  % Nombre de fonctions de base solides
  nbnnif = 3;  % Nombre de fonctions de base fluides

  if strcmp(support1,liste_intg{1})
%   RIGIDITE

%    nbptg  = 4; % Nombre de points d'integration
%%   Version poids central negatif Zienckiewick p. 201
%%   A ne pas utiliser en non-lineaire
%    disp('Modifier integration negative')
%    COOR   = [1./3. 1./3.
%              1./5. 1./5.
%              3./5. 1./5.
%              1./5. 3./5.];
%    WEIGHT = [-27./96.  25./96.  25./96. 25./96.];

%   Version 3 points de Gauss a l'interieur de l'element
%   Il existe aussi une version 3 points d'integration au milieu
%   des aretes (ordre a confirmer).
    nbptg  = 3; % cf Zienkiewicz
    [COOR,WEIGHT] = TRI_Integr(nbptg);
%%    COOR   = [1./6. 1./6.
%%              2./3. 1./6.
%%              1./6. 2./3.];
%%    WEIGHT = [1./6. 1./6. 1./6.];

  else
    support1
    error('support not implemented')
  end
elseif strcmp(mode1,'AXIS')
% AXIS Axisymmetric

  idimr  = 2;  % Nombre de coordonnees
  nbnnis = 6;  % Nombre de fonctions de base solides
  nbnnif = 3;  % Nombre de fonctions de base fluides

  if strcmp(support1,liste_intg{1})
%   RIGIDITE

    nbptg  = 7; % cf Castem et Bathe p. 280
    [COOR,WEIGHT] = TRI_Integr(nbptg);
%%    COOR   = [X101 X101
%%              X797 X101
%%              X101 X797
%%              X470 X059
%%              X470 X470
%%              X059 X470
%%              1./3. 1./3.];
%%    WEIGHT = [P125*0.5 P125*0.5 P125*0.5 P132*0.5 P132*0.5 P132*0.5 .1125D0];

  else
    support1
    error('support not implemented')
  end
else
  mode1
  error('mode not implemented')
end

% Values of shape functions and shape function derivatives
% at integration points
  PHI    = zeros(nbnnis+nbnnif,nbptg);
  DPHI   = zeros(idimr,nbnnis+nbnnif,nbptg);
  for ptg = 1:nbptg
    x = COOR(ptg,1); y = COOR(ptg,2);
    PHI(:,ptg) = [(1.-x-y)*(1.-2.*x-2.*y)
                  x*(2.*x-1.)
                  y*(2.*y-1.)
                  4.*x*(1.-x-y)
                  4.*x*y
                  4.*y*(1.-x-y)
                  1.-x-y
                  x
                  y];
    DPHI(:,:,ptg) = [(4.*x+4.*y-3)  (4.*x+4.*y-3)
                     (4.*x-1.)      0.
                     0.             (4.*y-1.)
                     (4.-8.*x-4.*y) (-4.*x)
                     (4.*y)         (4.*x)
                     (-4.*y)        (4.-4.*x-8.*y)
                     -1.            -1.
                     1.             0.
                     0.             1.]';
  end
