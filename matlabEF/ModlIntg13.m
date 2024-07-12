function [modl1,intg1] = ModlIntg13(mail1,mot1,support1,mode1, ...
                                    idim,varargin)
% Select model and integration information
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 01 / 08 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 02 / 09 / 2002
%  Ajout des deformations/contraintes
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 11 / 2002
%  Ajout du modele poreux, pas de duplication des fcts de base,
%  ajout TRI6 isoparametrique
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 11 / 2002
%  modification des modeles version 1.3 passe a version 1.4
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 22 / 11 / 2002
%  ajout du modele permeabilite
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 12 / 2002
%  ajout des elements SEG2, SEG3
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 13 / 12 / 2002
%  reunification de la formulation fluide
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 21 / 01 / 2003
%  ajout POI1
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 05 / 02 / 2003
%  ajout du mode d'analyse BARR (elements de structure)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 05 / 2003
%  possibilite de surcharger les points d'integration
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 23 / 09 / 2003
%  ajout de la formulation thermique, identique a fluide
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 12 / 2003
%  ajout du mode d'analyse TIMO (elements de structure poutre de
%  Timoshenko a courbure negligee)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 12 / 2003
%  ajout du mode d'analyse POUT (elements de structure poutre d'Euler
%  Bernoulli a courbure negligee)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 04 / 2004
%  ajout QUA4 fluide
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 26 / 04 / 2004
%  passage au 3D pour les barres et POI1
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 11 / 07 / 2004
%  ajout du mode d'analyse DKIR (elements de structure plaque de
%  Kirchhoff discret suppose plans)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 13 / 12 / 2004
%  ajout QUA8
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 05 / 2005
%  ajout du mode d'analyse JOIN (elements de joint)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 02 / 08 / 2005
%  ajout de la capacite thermique en QUA4
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 09 / 2005
%  ajout du QUA8
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 17 / 01 / 2006
%  passage au 3D pour les POUT
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 07 / 07 / 2006
%  Correction de bug pour BARR 3D
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 24 / 08 / 2006
%  Ajout TET4
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 31 / 08 / 2007
%  Ajout formulation TRANSPORT
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 29 / 08 / 2007
%  Ajout PRI6
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 08 / 2007
%  Ajout AXIS
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 24 / 11 / 2007
%  Ajout POREUX/TET5b, TE10
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 07 / 2008
%  Ajout ELASTIQUE, TRID, SEG2, TRI3 (composantes 3D, repere global)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 09 / 08 / 2008
%  Ajout POREUX/TET9b
% DUREISSEIX David  LaMCoS                           le 08 / 10 / 2010
%  ajout du mode d'analyse DSHE (elements de structure plaque avec
%  cisaillement discret suppose plans)
% DUREISSEIX David  LaMCoS                           le 13 / 06 / 2014
%  ajout THERMIQUE,BARR (thermique 1D)
% DUREISSEIX David  LaMCoS                           le 02 / 10 / 2019
%  Ajout ELASTIQUE, TRID, CU20
%
% Fournit le modele modl1 et le segment d'integration intg1
%
% [modl1,intg1] = ModlIntg13(mail1,mot1,support1,mode1,idim)
% [modl1,intg1] = ModlIntg13(mail1,mot1,support1,mode1,idim,intg2)
%
% Entrees
%   mail1	maillage
%   mot1	formulation (ELASTIQUE,POREUX,FLUIDE,THERMIQUE,
%                            TRANSPORT)
%   support1	support d'integration 
%               (RIGIDITE,MASSE,COMPRESSIBILITE,PERMEABILITE,
%                ADVECTION,FLUX)
%   mode1	mode d'analyse
%               (AXIS,COPL,DEPL,BARR,TIMO,POUT,DKIR,DSHE,JOIN)
%   idim	dimension de l'espace physique (2 ou 3)
%
% Entrees optionnelles
%   intg2	segment d'integration eventuel qui va etre utilise
%               pour imposer les points et poids d'integration
% Sorties
%   modl1	modele
%   intg1	segment d'integration
%
% Actuellement, il est possible de faire
%   ELASTIQUE
%     ELASTIQUE, TRID
%       ELASTIQUE, TRID, RIGIDITE
%       ELASTIQUE, TRID, MASSE
%       POI1 (ELASTIQUE, TRID)
%       SEG2 (ELASTIQUE, TRID)
%       TRI3 (ELASTIQUE, TRID)
%       CUB8 (ELASTIQUE, TRID)
%       QUA4 (ELASTIQUE, TRID)
%       PRI6 (ELASTIQUE, TRID)
%       TET4 (ELASTIQUE, TRID)
%       TE10 (ELASTIQUE, TRID)
%       CU20 (ELASTIQUE, TRID)
%     ELASTIQUE, (COPL Plane stress | DEPL Plane strain)
%       ELASTIQUE, (COPL | DEPL), RIGIDITE
%       ELASTIQUE, (COPL | DEPL), MASSE
%       QUA4 (ELASTIQUE, COPL | DEPL)
%       TRI3 (ELASTIQUE, COPL | DEPL)
%       TRI6 (ELASTIQUE, COPL | DEPL)
%       SEG2 (ELASTIQUE, COPL | DEPL)
%       SEG3 (ELASTIQUE, COPL | DEPL)
%       POI1 (ELASTIQUE, COPL | DEPL)
%       QUA8 (ELASTIQUE, COPL | DEPL)
%     ELASTIQUE, AXIS Axisymmetric
%       ELASTIQUE, AXIS, RIGIDITE
%       ELASTIQUE, AXIS, MASSE
%       SEG2 (ELASTIQUE, AXIS)
%       TRI3 (ELASTIQUE, AXIS)
%       TRI6 (ELASTIQUE, AXIS)
%     BARR Barre (element de structure avec base locale)
%       ELASTIQUE, BARR, RIGIDITE
%       ELASTIQUE, BARR, MASSE
%       SEG2 (ELASTIQUE, BARR)
%     TIMO Poutre de Timoshenko, courbure negligee
%       ELASTIQUE, TIMO, RIGIDITE
%       SEG2 (ELASTIQUE, TIMO) poutre droite de Timoshenko
%     POUT Poutre d'Euler Bernoulli, courbure negligee
%       ELASTIQUE, POUT, RIGIDITE
%       ELASTIQUE, POUT, MASSE
%       SEG2 (ELASTIQUE, POUT) poutre droite d'Euler-Bernoulli
%     DKIR Plaque de Kirchhoff discret, supposee plane
%       ELASTIQUE, DKIR, RIGIDITE
%       QUA4 (ELASTIQUE, DKIR) plaque plane de Kirchhoff discret
%       TRI3 (ELASTIQUE, DKIR) plaque plane de Kirchhoff discret
%     JOIN Modele de joint (epaisseur nulle)
%       ELASTIQUE, JOIN, RIGIDITE
%       RAC2 (ELASTIQUE, JOIN) Element de joint 2D
%       RAC3 (ELASTIQUE, JOIN) Element de joint 2D
%   POREUX
%     POREUX, (COPL | DEPL | AXIS)
%       POREUX, (COPL | DEPL | AXIS), RIGIDITE
%       TRI6 a bords droits (POREUX, COPL | DEPL | AXIS)
%     POREUX, TRID
%       TET4 Interpolation degree 1 for both displacement/pressure
%       TET5b Interpolation degree 1 + bubble (pyramid) for displacement
%       TET9b Interpolation degree 1 + bubble (pyramid) for displacement
%       TE10 Straight-edge tetrahedron with 10 nodes
%   FLUIDE
%     FLUIDE, (COPL | DEPL | AXIS)
%       FLUIDE, (COPL | DEPL | AXIS), COMPRESSIBILITE
%       FLUIDE, (COPL | DEPL | AXIS), PERMEABILITE
%       QUA4 (FLUIDE, COPL | DEPL)
%       TRI3 (FLUIDE, COPL | DEPL | AXIS)
%       TRI6 (FLUIDE, COPL | DEPL | AXIS)
%       SEG2 (FLUIDE, COPL | DEPL)
%     FLUIDE, DKIR
%       FLUIDE, DKIR, COMPRESSIBILITE
%       DKT (FLUIDE, DKIR)
%     FLUIDE, TRID
%       FLUIDE, TRID, COMPRESSIBILITE
%       FLUIDE, TRID, PERMEABILITE
%       SEG2 (FLUIDE, TRID)
%       TET4 (FLUIDE, TRID)
%       CUB8 (FLUIDE, TRID)
%       PRI6 (FLUIDE, TRID)
%     FLUIDE, BARR
%       FLUIDE, BARR, COMPRESSIBILITE
%       FLUIDE, BARR, PERMEABILITE
%       SEG2 (FLUIDE, BARR)
%   THERMIQUE (seuls les noms changent / FLUIDE)
%     THERMIQUE, (COPL | DEPL | AXIS)
%       THERMIQUE, (COPL | DEPL | AXIS), CAPACITE
%       THERMIQUE, (COPL | DEPL | AXIS), CONDUCTIVITE
%       TRI3 (THERMIQUE, COPL | DEPL)
%       TRI6 (THERMIQUE, COPL | DEPL)
%       SEG2 (THERMIQUE, COPL | DEPL)
%       QUA4 (THERMIQUE, COPL | DEPL)
%     THERMIQUE, TRID
%       THERMIQUE, TRID, CAPACITE
%       THERMIQUE, TRID, CONDUCTIVITE
%       TET4 (THERMIQUE, TRID)
%   TRANSPORT
%     TRANSPORT, (COPL | DEPL)
%       TRANSPORT, (COPL | DEPL), MASSE
%       TRANSPORT, (COPL | DEPL), ADVECTION
%       TRANSPORT, (COPL | DEPL), FLUX
%       TRI3 (TRANSPORT, COPL | DEPL)

% Extrait de la doc castem, pour info :
% Noms de degres de liberte
% Inconnue : UX UY UZ UT RX RY RZ RT RR P  PI  T TSUP TINF LX  TH   FC
% Duale    : FX FY FZ FT MX MY MZ MT MR FP FPI Q QSUP QINF FLX FLUX ED
% Inconnue : UR UZ P
% Duale    : FR FZ FP
% Noms de composantes
% Deformation : EPSS EPXX EPYY EPZZ GAXY DRSN DRN
% Contrainte  : EFFX SMXX SMYY SMZZ TAXY SMSN SMN
% Deformation : EPRR EPZZ EPTT GARZ
% Contrainte  : SMRR SMZZ SMTT TARZ

narg = nargin; % Nombre d'arguments
if (narg == 5)
  clear intg2;
elseif (narg == 6)
  intg2 = varargin{1};
else
  error('Bad number of arguments')
end

GlobalVar;

nbzone1 = length(mail1);
if exist('intg2')
  if (nbzone1 ~= length(intg2))
    nbzone1
    length(intg2)
    error('Bad number of subzones')
  end
end

clear ddlploc ddldloc; % Local basis dof, if there are any
clear DPHIN;           % Basis function derivates at nodes

clear modl1;
clear intg1;
for zo1 = 1:nbzone1
  type1 = mail1{zo1}.TYPE;

  if strcmp(mot1,liste_model{1})

%   ELASTIQUE
%   """""""""
    if strcmp(mode1,liste_mode{1})
%     ELASTIQUE, TRID
      if (idim ~= 3)
        mode1
        idim
        error('Only 3D for this mode')
      end
      ddlp = liste_ddlp;       % noms des ddl primaux
      ddld = liste_ddld;       % noms des ddl duaux
      if strcmp(support1,liste_intg{1})
%       ELASTIQUE, TRID, RIGIDITE
        comp     = liste_strain;     % noms des composantes primales
        comd     = liste_stress;     % noms des composantes duales
      elseif strcmp(support1,liste_intg{2})
%       ELASTIQUE, TRID, MASSE
        comp = liste_ddlp;     % noms des composantes primales
        comd = liste_ddld;     % noms des composantes duales
      else
        mot1
        mode1
        support1
        error('type of support not implemented in 3D')
      end
      if strcmp(type1,'POI1')
%       POI1 (ELASTIQUE, TRID)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = POI1_SegmentIntegr3(mode1,intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = POI1_SegmentIntegr3(mode1);
        end
        nddp = [1 2 3]; % ddl primaux
        nddd = nddp;    % ddl duaux
        nnop = [1 1 1]; % noeuds primaux
        nnod = nnop;   % noeuds duaux
        nnip = [];     % fct de base primales
        nnid = [];     % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [];    % noeuds pour transformation
        nnit1 = [];    % fct de base pour transformation
      elseif strcmp(type1,'SEG2')
%       SEG2 (ELASTIQUE, TRID)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 2 3 1 2 3];    % ddl primaux
        nddd = nddp;             % ddl duaux
        nnop = [1 1 1 2 2 2];    % noeuds primaux
        nnod = nnop;             % noeuds duaux
        nnip = nnop;             % fct de base primales
        nnid = nnod;             % fct de base duales
        ncop = [1:length(comp)]; % numeros des composantes primales
        ncod = [1:length(comd)]; % numeros des composantes duales
        nnot1 = [1 2];           % noeuds pour transformation
        nnit1 = [1 2];           % fct de base pour transformation
      elseif strcmp(type1,'TRI3')
%       TRI3 (ELASTIQUE, TRID)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 2 3 1 2 3 1 2 3];    % ddl primaux
        nddd = nddp;                   % ddl duaux
        nnop = [1 1 1 2 2 2 3 3 3];    % noeuds primaux
        nnod = nnop;                   % noeuds duaux
        nnip = nnop;                   % fct de base primales
        nnid = nnod;                   % fct de base duales
        ncop = [1:length(comp)]; % numeros des composantes primales
        ncod = [1:length(comd)]; % numeros des composantes duales
        nnot1 = [1 2 3];         % noeuds pour transformation
        nnit1 = [1 2 3];         % fct de base pour transformation
      elseif strcmp(type1,'CUB8')
%       CUB8 (ELASTIQUE, TRID)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = CUB8_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = CUB8_SegmentIntegr3(support1,mode1);
        end
%       ddl primaux
        nddp = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3];
%       ddl duaux
        nddd = nddp;
%       noeuds primaux
        nnop = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
%       noeuds duaux
        nnod = nnop;
%       fct de base primales
        nnip = nnop;
%       fct de base duales
        nnid = nnod;
%       numeros des composantes primales
        ncop = [1:length(comp)];
%       numeros des composantes duales
        ncod = [1:length(comd)];
%       noeuds pour transformation
        nnot1 = [1 2 3 4 5 6 7 8];
%       fct de base pour transformation
        nnit1 = [1 2 3 4 5 6 7 8];
      elseif strcmp(type1,QUA4)
%       QUA4 (ELASTIQUE, TRID)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = QUA4_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = QUA4_SegmentIntegr3(support1,mode1);
        end
%       ddl primaux
        nddp = [1 2 3 1 2 3 1 2 3 1 2 3];
%       ddl duaux
        nddd = nddp;
%       noeuds primaux
        nnop = [1 1 1 2 2 2 3 3 3 4 4 4];
%       noeuds duaux
        nnod = nnop;
%       fct de base primales
        nnip = nnop;
%       fct de base duales
        nnid = nnod;
%       numeros des composantes primales
        ncop = [1:length(comp)];
%       numeros des composantes duales
        ncod = [1:length(comd)];
%       noeuds pour transformation
        nnot1 = [1 2 3 4];
%       fct de base pour transformation
        nnit1 = [1 2 3 4];
      elseif strcmp(type1,'PRI6')
%       PRI6 (ELASTIQUE, TRID)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = PRI6_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = PRI6_SegmentIntegr3(support1,mode1);
        end
%       ddl primaux
        nddp = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3];
%       ddl duaux
        nddd = nddp;
%       noeuds primaux
        nnop = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6];
%       noeuds duaux
        nnod = nnop;
%       fct de base primales
        nnip = nnop;
%       fct de base duales
        nnid = nnod;
%       numeros des composantes primales
        ncop = [1:length(comp)];
%       numeros des composantes duales
        ncod = [1:length(comd)];
%       noeuds pour transformation
        nnot1 = [1 2 3 4 5 6];
%       fct de base pour transformation
        nnit1 = [1 2 3 4 5 6];
      elseif strcmp(type1,'TET4')
%       TET4 (ELASTIQUE, TRID)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TET4_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TET4_SegmentIntegr3(support1,mode1);
        end
%       ddl primaux
        nddp = [1 2 3 1 2 3 1 2 3 1 2 3];
%       ddl duaux
        nddd = nddp;
%       noeuds primaux
        nnop = [1 1 1 2 2 2 3 3 3 4 4 4];
%       noeuds duaux
        nnod = nnop;
%       fct de base primales
        nnip = nnop;
%       fct de base duales
        nnid = nnod;
%       numeros des composantes primales
        ncop = [1:length(comp)];
%       numeros des composantes duales
        ncod = [1:length(comd)];
%       noeuds pour transformation
        nnot1 = [1 2 3 4];
%       fct de base pour transformation
        nnit1 = [1 2 3 4];
      elseif strcmp(type1,'TE10')
%       TE10 (ELASTIQUE, TRID)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TE10_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TE10_SegmentIntegr3(support1,mode1);
        end
%       ddl primaux
        nddp = [1 2 3 1 2 3 1 2 3 1 2 3 ...
                1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3];
%       ddl duaux
        nddd = nddp;
%       noeuds primaux
        nnop = [1 1 1 2 2 2 3 3 3 4 4 4 ...
                5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10];
%       noeuds duaux
        nnod = nnop;
%       fct de base primales
        nnip = nnop;
%       fct de base duales
        nnid = nnod;
%       numeros des composantes primales
        ncop = [1:length(comp)];
%       numeros des composantes duales
        ncod = [1:length(comd)];
%       noeuds pour transformation
        nnot1 = [1 2 3 4 5 6 7 8 9 10];
%       fct de base pour transformation
        nnit1 = [1 2 3 4 5 6 7 8 9 10];
      elseif strcmp(type1,'CU20')
%       CU20 (ELASTIQUE, TRID)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = CU20_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = CU20_SegmentIntegr3(support1,mode1);
        end
%       ddl primaux
        nddp = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 ...
                1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3];
%       ddl duaux
        nddd = nddp;
%       noeuds primaux
        nnop = [1 1 1 2 2 2 3 3 3 4 4 4 ...
                5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 ...
                11 11 11 12 12 12 13 13 13 14 14 14 ...
                15 15 15 16 16 16 17 17 17 18 18 18 19 19 19 20 20 20];
%       noeuds duaux
        nnod = nnop;
%       fct de base primales
        nnip = nnop;
%       fct de base duales
        nnid = nnod;
%       numeros des composantes primales
        ncop = [1:length(comp)];
%       numeros des composantes duales
        ncod = [1:length(comd)];
%       noeuds pour transformation (isoparametrique)
        nnot1 = [1:20];
%       fct de base pour transformation
        nnit1 = [1:20];
      else
        mot1
        mode1
        type1
        error('type of element not implemented-0')
      end
    elseif (strcmp(mode1,liste_mode{2}) || strcmp(mode1,liste_mode{3}))
%     ELASTIQUE, (COPL Plane stress | DEPL Plane strain)
      if (idim ~= 2)
        mode1
        idim
        error('Only 2D for this mode')
      end
      ddlp = liste_ddlp2;       % noms des ddl primaux
      ddld = liste_ddld2;       % noms des ddl duaux
      if strcmp(support1,liste_intg{1})
%       ELASTIQUE, (COPL | DEPL), RIGIDITE
        comp = liste_strain2;     % noms des composantes primales
        comd = liste_stress2;     % noms des composantes duales
      elseif strcmp(support1,liste_intg{2})
%       ELASTIQUE, (COPL | DEPL), MASSE
        comp = liste_ddlp2;     % noms des composantes primales
        comd = liste_ddld2;     % noms des composantes duales
      else
        mot1
        mode1
        support1
        error('type of support not implemented')
      end
      if strcmp(type1,list_type_C3M1{1})
%       QUA4 (ELASTIQUE, COPL | DEPL)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = QUA4_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = QUA4_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 2 1 2 1 2 1 2]; % numeros des ddl primaux
        nddd = nddp;              % numeros des ddl duaux
        nnop = [1 1 2 2 3 3 4 4]; % noeuds primaux
        nnod = nnop;              % noeuds duaux
        nnip = nnop;              % fct de base primales
        nnid = nnod;              % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4];        % noeuds pour transformation
        nnit1 = [1 2 3 4];        % fct de base pour transformation
      elseif strcmp(type1,list_type_C3M1{2})
%       TRI3 (ELASTIQUE, COPL | DEPL)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 2 1 2 1 2]; % ddl primaux
        nddd = nddp;          % ddl duaux
        nnop = [1 1 2 2 3 3]; % noeuds primaux
        nnod = nnop;          % noeuds duaux
        nnip = nnop;          % fct de base primales
        nnid = nnod;          % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3];      % noeuds pour transformation
        nnit1 = [1 2 3];      % fct de base pour transformation
      elseif strcmp(type1,list_type_C3M2{2})
%       TRI6 (ELASTIQUE, COPL | DEPL)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TRI6_SegmentIntegr2(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TRI6_SegmentIntegr2(support1,mode1);
        end
        nddp = [1 2 1 2 1 2 1 2 1 2 1 2]; % ddl primaux
        nddd = nddp;                      % ddl duaux
        nnop = [1 1 2 2 3 3 4 4 5 5 6 6]; % noeuds primaux
        nnod = nnop;                      % noeuds duaux
        nnip = nnop;                      % fct de base primales
        nnid = nnod;                      % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4 5 6];    % noeuds pour transformation
        nnit1 = [1 2 3 4 5 6];    % fct de base pour transformation
      elseif strcmp(type1,list_type_C3M1{3})
%       SEG2 (ELASTIQUE, COPL | DEPL)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 2 1 2]; % ddl primaux
        nddd = nddp;      % ddl duaux
        nnop = [1 1 2 2]; % noeuds primaux
        nnod = nnop;      % noeuds duaux
        nnip = nnop;      % fct de base primales
        nnid = nnod;      % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2];            % noeuds pour transformation
        nnit1 = [1 2];            % fct de base pour transformation
      elseif strcmp(type1,list_type_C3M2{3})
%       SEG3 (ELASTIQUE, COPL | DEPL)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = SEG3_SegmentIntegr2(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = SEG3_SegmentIntegr2(support1,mode1);
        end
        nddp = [1 2 1 2 1 2]; % ddl primaux
        nddd = nddp;          % ddl duaux
        nnop = [1 1 2 2 3 3]; % noeuds primaux
        nnod = nnop;          % noeuds duaux
        nnip = nnop;          % fct de base primales
        nnid = nnod;          % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3];          % noeuds pour transformation
        nnit1 = [1 2 3];          % fct de base pour transformation
      elseif strcmp(type1,list_type_C3M1{4})
%       POI1 (ELASTIQUE, COPL | DEPL)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT]= POI1_SegmentIntegr3(mode1,intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = POI1_SegmentIntegr3(mode1);
        end
        nddp = [1 2]; % ddl primaux
        nddd = nddp;  % ddl duaux
        nnop = [1 1]; % noeuds primaux
        nnod = nnop;  % noeuds duaux
        nnip = [];    % fct de base primales
        nnid = [];    % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [];   % noeuds pour transformation
        nnit1 = [];   % fct de base pour transformation
      elseif strcmp(type1,'QUA8')
%       QUA8 (ELASTIQUE, COPL | DEPL)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = QUA8_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = QUA8_SegmentIntegr3(support1,mode1);
        end
	% numeros des ddl primaux
        nddp = [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2];
        nddd = nddp;              % numeros des ddl duaux
	% noeuds primaux
        nnop = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];
        nnod = nnop;              % noeuds duaux
        nnip = nnop;              % fct de base primales
        nnid = nnod;              % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4 5 6 7 8]; % noeuds pour transformation
        nnit1 = [1 2 3 4 5 6 7 8]; % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-1')
      end

    elseif strcmp(mode1,'AXIS')
%     ELASTIQUE, AXIS Axisymmetric
      if (idim ~= 2)
        mode1
        idim
        error('Only 2D for this mode')
      end
      ddlp = liste_ddlp2a;       % noms des ddl primaux
      ddld = liste_ddld2a;       % noms des ddl duaux
      if strcmp(support1,liste_intg{1})
%       ELASTIQUE, AXIS, RIGIDITE
        comp = liste_strain2a;     % noms des composantes primales
        comd = liste_stress2a;     % noms des composantes duales
      elseif strcmp(support1,liste_intg{2})
%       ELASTIQUE, AXIS, MASSE
        comp = liste_ddlp2a;     % noms des composantes primales
        comd = liste_ddld2a;     % noms des composantes duales
      else
        mot1
        mode1
        support1
        error('type of support not implemented')
      end
      if strcmp(type1,'SEG2')
%       SEG2 (ELASTIQUE, AXIS)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 2 1 2]; % ddl primaux
        nddd = nddp;      % ddl duaux
        nnop = [1 1 2 2]; % noeuds primaux
        nnod = nnop;      % noeuds duaux
        nnip = nnop;      % fct de base primales
        nnid = nnod;      % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2];            % noeuds pour transformation
        nnit1 = [1 2];            % fct de base pour transformation
      elseif strcmp(type1,'TRI3')
%       TRI3 (ELASTIQUE, AXIS)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 2 1 2 1 2]; % ddl primaux
        nddd = nddp;          % ddl duaux
        nnop = [1 1 2 2 3 3]; % noeuds primaux
        nnod = nnop;          % noeuds duaux
        nnip = nnop;          % fct de base primales
        nnid = nnod;          % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3];      % noeuds pour transformation
        nnit1 = [1 2 3];      % fct de base pour transformation
      elseif strcmp(type1,'TRI6')
%       TRI6 (ELASTIQUE, AXIS)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TRI6_SegmentIntegr2(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TRI6_SegmentIntegr2(support1,mode1);
        end
        nddp = [1 2 1 2 1 2 1 2 1 2 1 2]; % ddl primaux
        nddd = nddp;                      % ddl duaux
        nnop = [1 1 2 2 3 3 4 4 5 5 6 6]; % noeuds primaux
        nnod = nnop;                      % noeuds duaux
        nnip = nnop;                      % fct de base primales
        nnid = nnod;                      % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4 5 6];    % noeuds pour transformation
        nnit1 = [1 2 3 4 5 6];    % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-1a')
      end

    elseif strcmp(mode1,liste_mode{4})
%     BARR Barre (element de structure avec base locale)
      if strcmp(support1,liste_intg{1})
%       ELASTIQUE, BARR, RIGIDITE
	switch idim
	  case 2,
            ddlp = liste_ddlp2;       % noms des ddl primaux
            ddld = liste_ddld2;       % noms des ddl duaux
            ddlploc = [{'U1'}];       % noms des ddl primaux locaux
            ddldloc = [{'F1'}];       % noms des ddl duaux locaux
	  case 3,
            ddlp = liste_ddlp;        % noms des ddl primaux
            ddld = liste_ddld;        % noms des ddl duaux
            ddlploc = [{'U1'}];       % noms des ddl primaux locaux
            ddldloc = [{'F1'}];       % noms des ddl duaux locaux
        end
        comp = liste_barr_strain;     % noms des composantes primales
        comd = liste_barr_stress;     % noms des composantes duales
      elseif strcmp(support1,liste_intg{2})
%       ELASTIQUE, BARR, MASSE
	switch idim
	  case 2,
            ddlp = liste_ddlp2;           % noms des ddl primaux
            ddld = liste_ddld2;           % noms des ddl duaux
            ddlploc = liste_local_ddlp2;  % noms des ddl primaux locaux
            ddldloc = liste_local_ddld2;    % noms des ddl duaux locaux
            comp = liste_local_ddlp2;     % noms des composantes primales
            comd = liste_local_ddld2;     % noms des composantes duales
	  case 3,
            ddlp = liste_ddlp;                % noms des ddl primaux
            ddld = liste_ddld;                % noms des ddl duaux
            ddlploc = liste_local_ddlp; % noms des ddl primaux locaux
            ddldloc = liste_local_ddld; % noms des ddl duaux locaux
            comp = liste_local_ddlp;          % noms des composantes primales
            comd = liste_local_ddld;          % noms des composantes duales
	end
      else
        mot1
        mode1
        support1
        error('type of support not implemented')
      end
      if strcmp(type1,list_type_C3M1{3})
%       SEG2 (ELASTIQUE, BARR)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1);
        end
        [PHIN,DPHIN,COORN,WEIGHTN] =SEG2_SegmentIntegr3('NOEUDS',mode1);
        nddp = repmat([1:idim],1,2);          % ddl primaux base globale
        nddd = nddp;                          % ddl duaux base globale
        nnop = [ones(1,idim) 2*ones(1,idim)]; % noeuds prim. base globle
        nnod = nnop;                          % noeuds duaux base globle
	switch support1
	  case 'RIGIDITE',
            nddploc = [1 1];     % ddl primaux base locale
            ndddloc = nddploc;   % ddl duaux base locale
            nnoploc = [1 2];     % noeuds primaux base locale
            nnodloc = nnoploc;   % noeuds duaux base locale
            nnip = nnoploc;      % fct de base primales base locale
            nnid = nnodloc;      % fct de base duales base locale
	  case 'MASSE',
	    if (idim == 2)
              nddploc = [1 2 1 2]; % ddl primaux base locale
              ndddloc = nddploc;   % ddl duaux base locale
              nnoploc = [1 1 2 2]; % noeuds primaux base locale
              nnodloc = nnoploc;   % noeuds duaux base locale
              nnip = nnoploc;      % fct de base primales base locale
              nnid = nnoploc;      % fct de base duales base locale
	    elseif (idim == 3)
              nddploc = [1 2 3 1 2 3]; % ddl primaux base locale
              ndddloc = nddploc;       % ddl duaux base locale
              nnoploc = [1 1 1 2 2 2]; % noeuds primaux base locale
              nnodloc = nnoploc;       % noeuds duaux base locale
              nnip = nnoploc;          % fct de base primales base locale
              nnid = nnoploc;          % fct de base duales base locale
	    end
	end
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2];            % noeuds pour transformation
        nnit1 = [1 2];            % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-4')
      end

    elseif strcmp(mode1,liste_mode{5})
%     TIMO Poutre de Timoshenko, courbure negligee
%     (element de structure avec base locale)
      if (idim == 2)
        ddlp = liste_ddlpr2;           % noms des ddl primaux
        ddld = liste_ddldr2;           % noms des ddl duaux
        ddlploc = liste_local_ddlpr2;  % noms des ddl primaux locaux
        ddldloc = liste_local_ddldr2;  % noms des ddl duaux locaux
      elseif (idim == 3)
        ddlp = liste_ddlp2;           % noms des ddl primaux
        ddld = liste_ddld2;           % noms des ddl duaux
        ddlploc = liste_local_ddlp2;  % noms des ddl primaux locaux
        ddldloc = liste_local_ddld2;  % noms des ddl duaux locaux
        error('PAS ENCORE IMPLANTE 2, DESOLE')
      else
        mode1
        idim
        error('Bad dimension for this mode')
      end
      if strcmp(support1,liste_intg{1})
%       ELASTIQUE, TIMO, RIGIDITE 
        if (idim == 2)
          comp = list_pout_strain2;     % noms des composantes primales
          comd = list_pout_stress2;     % noms des composantes duales
        elseif (idim == 3)
          comp = list_pout_strain;     % noms des composantes primales
          comd = list_pout_stress;     % noms des composantes duales
        else
          mode1
          idim
          support1
          error('Bad dimension for this support and this mode')
        end
      else
        mot1
        mode1
        support1
        error('type of support not implemented')
      end
      if strcmp(type1,list_type_C3M1{3})
%       SEG2 (ELASTIQUE, TIMO) poutre droite de Timoshenko
%                              courbure eventuelle negligee !
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1);
        end
        [PHIN,DPHIN,COORN,WEIGHTN] =SEG2_SegmentIntegr3('NOEUDS',mode1);
        if (idim == 2)
          nddp = [1 2 3 1 2 3]; % ddl primaux base globale
          nddd = nddp;          % ddl duaux base globale
          nnop = [1 1 1 2 2 2]; % noeuds prim. base globle
          nnod = nnop;          % noeuds duaux base globle
        elseif (idim == 3)
          nddp = [1 2 3 4 5 6 1 2 3 4 5 6]; % ddl primaux base globale
          nddd = nddp;                      % ddl duaux base globale
          nnop = [1 1 1 1 1 1 2 2 2 2 2 2]; % noeuds prim. base globle
          nnod = nnop;                      % noeuds duaux base globle
        end
        nddploc = nddp;      % ddl primaux base locale
        ndddloc = nddploc;   % ddl duaux base locale
        nnoploc = nnop;      % noeuds primaux base locale
        nnodloc = nnoploc;   % noeuds duaux base locale
        nnip = nnoploc;      % fct de base primales base locale
        nnid = nnodloc;      % fct de base duales base locale
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2];            % noeuds pour transformation
        nnit1 = [1 2];            % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-4')
      end

    elseif strcmp(mode1,liste_mode{6})
%     POUT Poutre d'Euler Bernoulli, courbure negligee
%     (element de structure avec base locale)
      if (idim == 2)
        ddlp = liste_ddlpr2;           % noms des ddl primaux
        ddld = liste_ddldr2;           % noms des ddl duaux
        ddlploc = liste_local_ddlpr2;  % noms des ddl primaux locaux
        ddldloc = liste_local_ddldr2;  % noms des ddl duaux locaux
      elseif (idim == 3)
        ddlp = liste_ddlpr;           % noms des ddl primaux
        ddld = liste_ddldr;           % noms des ddl duaux
        ddlploc = liste_local_ddlpr;  % noms des ddl primaux locaux
        ddldloc = liste_local_ddldr;  % noms des ddl duaux locaux
%%DD        error('PAS ENCORE IMPLANTE 3D, DESOLE')
      else
        mode1
        idim
        error('Bad dimension for this mode')
      end
      if strcmp(support1,liste_intg{1})
%       ELASTIQUE, POUT, RIGIDITE
        if (idim == 2)
          comp = [{'EPS1'} {'C3'}];     % noms des composantes primales
          comd = [{'EFF1'} {'MOM3'}];   % noms des composantes duales
        elseif (idim == 3)
          comp = [{'EPS1'} {'C1'}   {'C2'}   {'C3'}];
          comd = [{'EFF1'} {'MOM1'} {'MOM2'} {'MOM3'}];
        else
          mode1
          idim
          support1
          error('Bad dimension for this support and this mode')
        end
      elseif strcmp(support1,'MASSE')
%       ELASTIQUE, POUT, MASSE
        if (idim == 2)
          comp = liste_local_ddlpr2;     % noms des composantes primales
          comd = liste_local_ddldr2;     % noms des composantes duales
        elseif (idim == 3)
          comp = liste_local_ddlpr;      % noms des composantes primales
          comd = liste_local_ddldr;      % noms des composantes duales
        end
      else
        mot1
        mode1
        support1
        error('type of support not implemented')
      end
      if strcmp(type1,list_type_C3M1{3})
%       SEG2 (ELASTIQUE, POUT) poutre droite d'Euler-Bernoulli
%                              courbure eventuelle negligee !
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT,DDPHI] = ...
             SEG2_SegmentIntegr3(support1,mode1,intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT,DDPHI] = ...
             SEG2_SegmentIntegr3(support1,mode1);
        end
        [PHIN,DPHIN,COORN,WEIGHTN] =SEG2_SegmentIntegr3('NOEUDS',mode1);
        if (idim == 2)
          nddp = [1 2 3 1 2 3];    % ddl primaux base globale
          nddd = nddp;             % ddl duaux base globale
          nnop = [1 1 1 2 2 2];    % noeuds prim. base globale
          nnod = nnop;             % noeuds duaux base globale
          nddploc = [1 2 3 1 2 3]; % ddl primaux base locale
          ndddloc = nddploc;       % ddl duaux base locale
          nnoploc = [1 1 1 2 2 2]; % noeuds primaux base locale
          nnodloc = nnoploc;       % noeuds duaux base locale
          nnip    = [5 1 3 6 2 4]; % fct de base primales base locale
          nnid    = nnip;          % fct de base duales base locale
        elseif (idim == 3)
          nddp = [1 2 3 4 5 6 1 2 3 4 5 6]; % ddl primaux base globale
          nddd = nddp;                      % ddl duaux base globale
          nnop = [1 1 1 1 1 1 2 2 2 2 2 2]; % noeuds prim. base globle
          nnod = nnop;                      % noeuds duaux base globle
%%DD          error('A faire en 3D...')
          nddploc = [1 2 3 4 5 6 1 2 3 4 5 6]; % ddl primaux base locale
          ndddloc = nddploc;                   % ddl duaux base locale
          nnoploc = [1 1 1 1 1 1 2 2 2 2 2 2]; % noeuds primaux base locale
          nnodloc = nnoploc;                   % noeuds duaux base locale
          nnip    = [5 1 1 5 3 3 6 2 2 6 4 4]; % fct de base prim. base locale
          nnid    = nnip;                      % fct de base duales base locale
        end
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2];            % noeuds pour transformation
        nnit1 = [5 6];            % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-4')
      end

    elseif strcmp(mode1,liste_mode{7})
%     DKIR Plaque de Kirchhoff discret, supposee plane
%     (element de structure avec base locale)
      if (idim == 2)
        error('DKIR PAS IMPLANTE EN 2D 2a, DESOLE')
      elseif (idim == 3)
%       noms des ddl primaux
        ddlp = [{'UX'} {'UY'} {'UZ'} {'RX'} {'RY'} {'RZ'}];
%       noms des ddl duaux correspondant
        ddld = [{'FX'} {'FY'} {'FZ'} {'MX'} {'MY'} {'MZ'}];
%       noms des ddl primaux locaux
        ddlploc = [{'U1'} {'U2'} {'U3'} {'R1'} {'R2'}];
%       noms des ddl duaux locaux correspondant
        ddldloc = [{'F1'} {'F2'} {'F3'} {'M1'} {'M2'}];
      else
        mode1
        idim
        error('Bad dimension for this mode')
      end
      if strcmp(support1,liste_intg{1})
%       ELASTIQUE, DKIR, RIGIDITE
        if (idim == 2)
          error('DKIR PAS IMPLANTE EN 2D 2a, DESOLE')
        elseif (idim == 3)
%         noms des composantes primales 
          comp = [{'EP11'} {'EP22'} {'GA12'} {'CP11'} {'CP22'} {'CG12'}];
%         noms des composantes duales correspondant 
          comd = [{'EF11'} {'EF22'} {'EG12'} {'MF11'} {'MF22'} {'MG12'}];
        else
          mode1
          idim
          support1
          error('Bad dimension for this support and this mode')
        end
      else
        mot1
        mode1
        support1
        error('type of support not implemented')
      end

      if strcmp(type1,list_type_C3M1{1})
%       QUA4 (ELASTIQUE, DKIR) plaque plane de Kirchhoff discret
%                              (DKQ for bending, QUA4 for membrane)
%                              courbure eventuelle negligee !
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT,DDPHI] = ...
             QUA4_SegmentIntegr3(support1,mode1,intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT,DDPHI] = ...
             QUA4_SegmentIntegr3(support1,mode1);
%          [PHI,DPHI,COOR,WEIGHT,DDPHI] = ...
%             QUA4_SegmentIntegr3('RIGIDITE12',mode1);
        end
        [PHIN,DPHIN,COORN,WEIGHTN] = QUA4_SegmentIntegr3('NOEUDS',mode1);
        if (idim == 3)
%         numeros des ddl primaux base globale
          nddp = [1 2 3 4 5 6]; % for first node
          nddp = [nddp nddp nddp nddp]; % for all nodes
%         numeros des ddl duaux base globale
          nddd = nddp;
%         numeros locaux des noeuds primaux base globale
          nnop = [1 1 1 1 1 1 2 2 2 2 2 2 ...
                  3 3 3 3 3 3 4 4 4 4 4 4];
%         numeros locaux des noeuds duaux base globale
          nnod = nnop;

%         numeros des ddl primaux base locale
          nddploc = [1 2 3 4 5]; % for first node
          nddploc = [nddploc nddploc nddploc nddploc]; % for all corner nodes
          nddploc = [nddploc 4 5 4 5 4 5 4 5]; % add ons for virtual edge nodes
%         numeros des ddl duaux base locale
          ndddloc = nddploc;
%         numeros locaux des noeuds primaux base locale
          nnoploc = [1 1 1 1 1 2 2 2 2 2 ...
                     3 3 3 3 3 4 4 4 4 4];
%         add ons for virtual edge nodes
          nnoploc = [nnoploc -1 -1 -2 -2 -3 -3 -4 -4];
%         numeros locaux des noeuds duaux base locale
          nnodloc = nnoploc;
%         numeros des fct de base primales base locale
          nnip    = [1 1 1 1 1 2 2 2 2 2 ...
                     3 3 3 3 3 4 4 4 4 4]; % for all corner nodes
          nnip    = [nnip 5 5 6 6 7 7 8 8]; % add ons for virtual edge nodes
%         numeros des fct de base duales base locale
          nnid    = nnip;
        end
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4];        % noeuds pour transformation
        nnit1 = [1 2 3 4];        % fct de base pour transformation

      elseif strcmp(type1,list_type_C3M1{2})
%       TRI3 (ELASTIQUE, DKIR) plaque plane de Kirchhoff discret
%                              (DKT for bending, TRI3 for membrane)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT,DDPHI] = ...
             TRI3_SegmentIntegr3(support1,mode1,intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT,DDPHI] = ...
             TRI3_SegmentIntegr3(support1,mode1);
        end
        [PHIN,DPHIN,COORN,WEIGHTN] = TRI3_SegmentIntegr3('NOEUDS',mode1);
        if (idim == 3)
%         numeros des ddl primaux base globale
          nddp = [1 2 3 4 5 6]; % for first node
          nddp = [nddp nddp nddp]; % for all nodes
%         numeros des ddl duaux base globale
          nddd = nddp;
%         numeros locaux des noeuds primaux base globale
          nnop = [1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3];
%         numeros locaux des noeuds duaux base globale
          nnod = nnop;

%         numeros des ddl primaux base locale
          nddploc = [1 2 3 4 5]; % for first node
          nddploc = [nddploc nddploc nddploc]; % for all corner nodes
          nddploc = [nddploc 4 5 4 5 4 5]; % add ons for virtual edge nodes
%         numeros des ddl duaux base locale
          ndddloc = nddploc;
%         numeros locaux des noeuds primaux base locale
          nnoploc = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3];
%         add ons for virtual edge nodes
          nnoploc = [nnoploc -1 -1 -2 -2 -3 -3];
%         numeros locaux des noeuds duaux base locale
          nnodloc = nnoploc;
%         numeros des fct de base primales base locale
          nnip    = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3]; % for all corner nodes
          nnip    = [nnip 4 4 5 5 6 6]; % add ons for virtual edge nodes
%         numeros des fct de base duales base locale
          nnid    = nnip;
        end
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3];        % noeuds pour transformation
        nnit1 = [1 2 3];        % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-4a')
      end

    elseif strcmp(mode1,liste_mode{9})
%     DSHE Plaque avec cisaillement discret, supposee plane
%     (element de structure avec base locale)
      if (idim == 3)
%       noms des ddl primaux
        ddlp = [{'UX'} {'UY'} {'UZ'} {'RX'} {'RY'} {'RZ'}];
%       noms des ddl duaux correspondant
        ddld = [{'FX'} {'FY'} {'FZ'} {'MX'} {'MY'} {'MZ'}];
%       noms des ddl primaux locaux
        ddlploc = [{'U1'} {'U2'} {'U3'} {'R1'} {'R2'}];
%       noms des ddl duaux locaux correspondant
        ddldloc = [{'F1'} {'F2'} {'F3'} {'M1'} {'M2'}];
      else
        mode1
        idim
        error('Bad dimension for this mode')
      end
      if strcmp(support1,liste_intg{1})
%       ELASTIQUE, DSHE, RIGIDITE
        if (idim == 3)
%         noms des composantes primales 
          comp = [{'EP11'} {'EP22'} {'GA12'} ...
                  {'CP11'} {'CP22'} {'CG12'} ...
                  {'G1'} {'G2'}];
%         noms des composantes duales correspondant 
          comd = [{'EF11'} {'EF22'} {'EG12'} ...
                  {'MF11'} {'MF22'} {'MG12'} ...
                  {'T1'} {'T2'}];
        else
          mode1
          idim
          support1
          error('Bad dimension for this support and this mode')
        end
      else
        mot1
        mode1
        support1
        error('type of support not implemented')
      end

      if strcmp(type1,list_type_C3M1{2})
%       TRI3 (ELASTIQUE, DSHE) plaque plane avec cisaillement discret
%                              (DST for bending, TRI3 for membrane)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT,DDPHI] = ...
             TRI3_SegmentIntegr3(support1,mode1,intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT,DDPHI] = ...
             TRI3_SegmentIntegr3(support1,mode1);
        end
        [PHIN,DPHIN,COORN,WEIGHTN] = TRI3_SegmentIntegr3('NOEUDS',mode1);
        if (idim == 3)
%         numeros des ddl primaux base globale
          nddp = [1 2 3 4 5 6]; % for first node
          nddp = [nddp nddp nddp]; % for all nodes
%         numeros des ddl duaux base globale
          nddd = nddp;
%         numeros locaux des noeuds primaux base globale
          nnop = [1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3];
%         numeros locaux des noeuds duaux base globale
          nnod = nnop;

%         numeros des ddl primaux base locale
          nddploc = [1 2 3 4 5]; % for first node
          nddploc = [nddploc nddploc nddploc]; % for all corner nodes
          nddploc = [nddploc 4 5 4 5 4 5]; % add ons for virtual edge nodes
%         numeros des ddl duaux base locale
          ndddloc = nddploc;
%         numeros locaux des noeuds primaux base locale
          nnoploc = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3];
%         add ons for virtual edge nodes (negative values)
          nnoploc = [nnoploc -1 -1 -2 -2 -3 -3];
%         numeros locaux des noeuds duaux base locale
          nnodloc = nnoploc;
%         numeros des fct de base primales base locale
          nnip    = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3]; % for all corner nodes
          nnip    = [nnip 4 4 5 5 6 6]; % add ons for virtual edge nodes
%         numeros des fct de base duales base locale
          nnid    = nnip;
        end
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3];        % noeuds pour transformation
        nnit1 = [1 2 3];        % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-4b')
      end

    elseif strcmp(mode1,liste_mode{8})
%     JOIN Modele de joint (epaisseur nulle)
%     (element avec base locale)
      if (idim == 2)
        ddlp = [{'UX'} {'UY'}];       % noms des ddl primaux
        ddld = [{'FX'} {'FY'}];       % noms des ddl duaux
        ddlploc = [{'U1'} {'U2'}];    % noms des ddl primaux locaux
        ddldloc = [{'F1'} {'F2'}];    % noms des ddl duaux locaux
      elseif (idim == 3)
        ddlp = [{'UX'} {'UY'} {'UZ'}];  % noms des ddl primaux
        ddld = [{'FX'} {'FY'} {'FZ'}];  % noms des ddl duaux
        ddlploc = [{'U1'} {'U2'} {'U3'}];  % noms des ddl primaux locaux
        ddldloc = [{'F1'} {'F2'} {'F3'}];  % noms des ddl duaux locaux
        error('PAS ENCORE IMPLANTE EN DIM 3, DESOLE')
      else
        mode1
        idim
        error('Bad dimension for this mode')
      end
      if strcmp(support1,liste_intg{1})
%       ELASTIQUE, JOIN, RIGIDITE
        if (idim == 2)
          comp = [{'DRSN'} {'DRN'}];   % noms des composantes primales
          comd = [{'SMSN'} {'SMN'}];   % noms des composantes duales
        else
          mode1
          idim
          support1
          error('Bad dimension for this support and this mode')
        end
      else
        mot1
        mode1
        support1
        error('type of support not implemented')
      end
      if strcmp(type1,'RAC2')
%       RAC2 (ELASTIQUE, JOIN) Element de joint 2D
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = ...
             RAC2_SegmentIntegr3(support1,mode1,intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = ...
             RAC2_SegmentIntegr3(support1,mode1);
        end
        [PHIN,DPHIN,COORN,WEIGHTN] = RAC2_SegmentIntegr3('NOEUDS',mode1);
        if (idim == 2)
%         numeros des ddl primaux base globale
          nddp = [1 2 1 2 1 2 1 2]; % numeros des ddl primaux
          nddd = nddp;              % numeros des ddl duaux
          nnop = [1 1 2 2 3 3 4 4]; % noeuds primaux
          nnod = nnop;              % noeuds duaux
          nnip = [1 1 2 2 2 2 1 1]; % fct de base primales
          nnid = nnip;              % fct de base duales
          ncop = [1:length(comp)];  % numeros des composantes primales
          ncod = [1:length(comd)];  % numeros des composantes duales
%          nnot1 = [1 2 3 4];        % noeuds pour transformation
%          nnit1 = [1 2 3 4];        % fct de base pour transformation
          nnot1 = [1 2];        % noeuds pour transformation
          nnit1 = [1 2];        % fct de base pour transformation
%         numeros des ddl primaux base locale
          nddploc = [1 2 1 2 1 2 1 2]; % ddl primaux base locale
          ndddloc = nddploc;           % ddl duaux base locale
          nnoploc = [1 1 2 2 3 3 4 4]; % noeuds primaux base locale
          nnodloc = nnoploc;           % noeuds duaux base locale
          nnip    = [1 1 2 2 2 2 1 1]; % fct de base prim. base locale
          nnid    = nnip;              % fct de base duales base locale
        end
      elseif strcmp(type1,'RAC3')
%       RAC3 (ELASTIQUE, JOIN) Element de joint 2D
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = RAC3_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = RAC3_SegmentIntegr3(support1,mode1);
        end
        [PHIN,DPHIN,COORN,WEIGHTN] = RAC3_SegmentIntegr3('NOEUDS',mode1);
        if (idim == 2)
%         numeros des ddl primaux base globale
          nddp = [1 2 1 2 1 2 1 2 1 2 1 2]; % numeros des ddl primaux
          nddd = nddp;                      % numeros des ddl duaux
          nnop = [1 1 2 2 3 3 4 4 5 5 6 6]; % noeuds primaux
          nnod = nnop;                      % noeuds duaux
          nnip = [1 1 2 2 3 3 3 3 2 2 1 1]; % fct de base primales
          nnid = nnip;                      % fct de base duales
          ncop = [1:length(comp)];  % numeros des composantes primales
          ncod = [1:length(comd)];  % numeros des composantes duales
%          nnot1 = [1 2 3 4];        % noeuds pour transformation
%          nnit1 = [1 2 3 4];        % fct de base pour transformation
          nnot1 = [1 2 3];        % noeuds pour transformation
          nnit1 = [1 2 3];        % fct de base pour transformation
%         numeros des ddl primaux base locale
          nddploc = [1 2 1 2 1 2 1 2 1 2 1 2]; % ddl primaux base locale
          ndddloc = nddploc;                   % ddl duaux base locale
          nnoploc = [1 1 2 2 3 3 4 4 5 5 6 6]; % noeuds primaux base locale
          nnodloc = nnoploc;                   % noeuds duaux base locale
          nnip    = [1 1 2 2 3 3 3 3 2 2 1 1]; % fct de base prim. base locale
          nnid    = nnip;                      % fct de base duales base locale
        end
      else
        mot1
        mode1
        type1
        error('type of element not implemented-6')
      end

    else
      mot1
      mode1
      error('mode not implemented 1')
    end

  elseif strcmp(mot1,liste_model{2})

%   POREUX
%   """"""
    if (strcmp(mode1,'COPL') || strcmp(mode1,'DEPL') || strcmp(mode1,'AXIS'))
%     POREUX, (COPL | DEPL | AXIS)
      if (idim ~= 2)
        mode1
        idim
        error('Only 2D for this mode... bis')
      end
      if (strcmp(mode1,'COPL') || strcmp(mode1,'DEPL'))
        ddlp = [liste_ddlp2 {'P'}];       % noms des ddl primaux
        ddld = [liste_ddld2 {'FP'}];      % noms des ddl duaux
        if strcmp(support1,liste_intg{1})
%         POREUX, (COPL | DEPL | AXIS), RIGIDITE
          comp = [liste_strain2 {'P'}];   % noms des composantes primales
          comd = [liste_stress2 {'FP'}];  % noms des composantes duales
        else
          mot1
          mode1
          support1
          error('type of support not implemented')
        end
      elseif strcmp(mode1,'AXIS')
        ddlp = [liste_ddlp2a {'P'}];       % noms des ddl primaux
        ddld = [liste_ddld2a {'FP'}];      % noms des ddl duaux
        if strcmp(support1,liste_intg{1})
%         POREUX, (COPL | DEPL | AXIS), RIGIDITE
          comp = [liste_strain2a {'P'}];   % noms des composantes primales
          comd = [liste_stress2a {'FP'}];  % noms des composantes duales
        else
          mot1
          mode1
          support1
          error('type of support not implemented')
        end
      end
      if strcmp(type1,list_type_C3M2{2})
%       TRI6 a bords droits (POREUX, COPL | DEPL | AXIS)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = ...
             TRI6d_Poreux_SegmentIntegr2(support1,mode1,intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = ...
             TRI6d_Poreux_SegmentIntegr2(support1,mode1);
        end
        nddp = [1 2 1 2 1 2 1 2 1 2 1 2 3 3 3]; % ddl primaux
        nddd = nddp;                            % ddl duaux
        nnop = [1 1 2 2 3 3 4 4 5 5 6 6 1 2 3]; % noeuds primaux
        nnod = nnop;                            % noeuds duaux
        nnip = [1 1 2 2 3 3 4 4 5 5 6 6 7 8 9]; % fct de base primales
        nnid = nnip;                            % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3];      % noeuds pour transformation
        nnit1 = [7 8 9];      % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-2')
      end
    elseif strcmp(mode1,'TRID')
%     POREUX, TRID
      if (idim ~= 3)
        mode1
        idim
        error('Only 3D for this mode...')
      end
      ddlp = [liste_ddlp {'P'}];       % noms des ddl primaux
      ddld = [liste_ddld {'FP'}];      % noms des ddl duaux
      if strcmp(support1,'RIGIDITE')
        comp = [liste_strain {'P'}];   % noms des composantes primales
        comd = [liste_stress {'FP'}];  % noms des composantes duales
      else
        mot1
        mode1
        support1
        error('type of support not implemented (TRID)')
      end

      if strcmp(type1,'TET4')
%       Interpolation degree 1 for both displacement/pressure
        disp('ModlIntg13: POREUX/TRID/TET4 en test - instable !!!!!!')
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = ...
             TET4_Poreux_SegmentIntegr3(support1,mode1,intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = ...
             TET4_Poreux_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 2 3 1 2 3 1 2 3 1 2 3 4 4 4 4]; % ddl primaux
        nddd = nddp;                              % ddl duaux
        nnop = [1 1 1 2 2 2 3 3 3 4 4 4 1 2 3 4]; % noeuds primaux
        nnod = nnop;                              % noeuds duaux
        nnip = [1 1 1 2 2 2 3 3 3 4 4 4 1 2 3 4]; % fct de base primales
        nnid = nnip;                              % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4];      % noeuds pour transformation
        nnit1 = [1 2 3 4];      % fct de base pour transformation
      elseif strcmp(type1,'TET5b')
%       Interpolation degree 1 + bubble (pyramid) for displacement
%       Classical interpolations changed for displacement
%       Simple degree 1 for pressure
%       Classical interpolations for pressure (and transformation)
%       Continuous interpolations (Taylor-Hood element)
        disp('ModlIntg13: POREUX/TRID/TET5b en test !!!!!!')
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = ...
             TET5b_Poreux_SegmentIntegr3(support1,mode1,intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = ...
             TET5b_Poreux_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 4 4 4 4]; % ddl primaux
        nddd = nddp;                                    % ddl duaux
        nnop = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 1 2 3 4]; % noeuds primaux
        nnod = nnop;                                    % noeuds duaux
        nnip = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 7 8 9]; % fct de base primales
        nnid = nnip;                                    % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4];      % noeuds pour transformation
        nnit1 = [6 7 8 9];      % fct de base pour transformation
      elseif strcmp(type1,'TET9b')
%       Interpolation degree 1 + bubble (pyramid) for displacement
%       Classical interpolations changed for displacement
%       Simple degree 1 for pressure
%       Classical interpolations for pressure (and transformation)
%       Continuous interpolations (Taylor-Hood element)
        disp('ModlIntg13: POREUX/TRID/TET9b en test !!!!!!')
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = ...
             TET9b_Poreux_SegmentIntegr3(support1,mode1,intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = ...
             TET9b_Poreux_SegmentIntegr3(support1,mode1);
        end
%       ddl primaux
        nddp = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 4 4 4 4];
        nddd = nddp;              % ddl duaux
%       noeuds primaux
        nnop = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 1 2 3 4];
        nnod = nnop;              % noeuds duaux
%       fct de base primales
        nnip = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 11 12 13];
        nnid = nnip;              % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4];        % noeuds pour transformation
        nnit1 = [10 11 12 13];    % fct de base pour transformation
      elseif strcmp(type1,'TE10')
%       Straight-edge tetrahedron with 10 nodes
%       Interpolation degree 2 for displacement, degree 1 for pressure
        disp('ModlIntg13: POREUX/TRID/TE10 en test !')
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = ...
             TE10_Poreux_SegmentIntegr3(support1,mode1,intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = ...
             TE10_Poreux_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 2 3 1 2 3 1 2 3 1 2 3 ...
                1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 ...
                4 4 4 4]; % ddl primaux
        nddd = nddp;      % ddl duaux
        nnop = [1 1 1 2 2 2 3 3 3 4 4 4 ...
                5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 ...
                1 2 3 4]; % noeuds primaux
        nnod = nnop;      % noeuds duaux
        nnip = [1 1 1 2 2 2 3 3 3 4 4 4 ...
                5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 ...
                11 12 13 14]; % fct de base primales
        nnid = nnip;          % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [ 1  2  3  4];    % noeuds pour transformation
        nnit1 = [11 12 13 14];    % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-POREUX/TRID')
      end
    else
      mot1
      mode1
      error('mode not implemented 2')
    end

  elseif strcmp(mot1,liste_model{3})

%   FLUIDE
%   """"""
    if (strcmp(mode1,'COPL') || strcmp(mode1,'DEPL') || strcmp(mode1,'AXIS'))
%     FLUIDE, (COPL | DEPL | AXIS)
      if (idim ~= 2)
        mode1
        idim
        error('Only 2D for this mode... ter')
      end
      ddlp = [{'P'}];       % noms des ddl primaux
      ddld = [{'FP'}];      % noms des ddl duaux
      if strcmp(support1,liste_intg{3})
%       FLUIDE, (COPL | DEPL | AXIS), COMPRESSIBILITE
        comp = [{'P'}];   % noms des composantes primales
        comd = [{'FP'}];  % noms des composantes duales
      elseif strcmp(support1,liste_intg{4})
%       FLUIDE, (COPL | DEPL | AXIS), PERMEABILITE
        if (strcmp(mode1,'COPL') || strcmp(mode1,'DEPL'))
          comp = [{'PX'} {'PY'}];     % noms des composantes primales
          comd = [{'WX'} {'WY'}];     % noms des composantes duales
        elseif strcmp(mode1,'AXIS')
          comp = [{'PR'} {'PZ'}];     % noms des composantes primales
          comd = [{'WR'} {'WZ'}];     % noms des composantes duales
        else
          mode1
          error('Mode not yet implemented... sorry')
        end
      else
        mot1
        mode1
        support1
        error('type of support not implemented')
      end
      if strcmp(type1,list_type_C3M1{1})
%       QUA4 (FLUIDE, COPL | DEPL)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = QUA4_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = QUA4_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1 1 1];         % ddl primaux
        nddd = nddp;              % ddl duaux
        nnop = [1 2 3 4];         % noeuds primaux
        nnod = nnop;              % noeuds duaux
        nnip = [1 2 3 4];         % fct de base primales
        nnid = nnip;              % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4];        % noeuds pour transformation
        nnit1 = [1 2 3 4];        % fct de base pour transformation
      elseif strcmp(type1,list_type_C3M1{2})
%       TRI3 (FLUIDE, COPL | DEPL | AXIS)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1 1]; % ddl primaux
        nddd = nddp;    % ddl duaux
        nnop = [1 2 3]; % noeuds primaux
        nnod = nnop;    % noeuds duaux
        nnip = [1 2 3]; % fct de base primales
        nnid = nnip;    % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3];      % noeuds pour transformation
        nnit1 = [1 2 3];      % fct de base pour transformation
      elseif strcmp(type1,list_type_C3M2{2})
%       TRI6 (FLUIDE, COPL | DEPL | AXIS)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TRI6_SegmentIntegr2(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TRI6_SegmentIntegr2(support1,mode1);
        end
        nddp = [1 1 1 1 1 1];     % ddl primaux
        nddd = nddp;              % ddl duaux
        nnop = [1 2 3 4 5 6];     % noeuds primaux
        nnod = nnop;              % noeuds duaux
        nnip = [1 2 3 4 5 6];     % fct de base primales
        nnid = nnod;              % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4 5 6];    % noeuds pour transformation
        nnit1 = [1 2 3 4 5 6];    % fct de base pour transformation
      elseif strcmp(type1,list_type_C3M1{3})
%       SEG2 (FLUIDE, COPL | DEPL)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1]; % ddl primaux
        nddd = nddp;  % ddl duaux
        nnop = [1 2]; % noeuds primaux
        nnod = nnop;  % noeuds duaux
        nnip = nnop;  % fct de base primales
        nnid = nnod;  % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2]; % noeuds pour transformation
        nnit1 = [1 2]; % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-3')
      end
    elseif strcmp(mode1,'DKIR')
%     FLUIDE, DKIR
      if (idim ~= 3)
        mode1
        idim
        error('Only 3D for this mode...')
      end
      ddlp = [{'P'}];       % noms des ddl primaux
      ddld = [{'FP'}];      % noms des ddl duaux
      if strcmp(support1,liste_intg{3})
%       FLUIDE, DKIR, COMPRESSIBILITE
        comp = [{'P'}];   % noms des composantes primales
        comd = [{'FP'}];  % noms des composantes duales
      else
        mot1
        mode1
        support1
        error('type of support not implemented 2aa')
      end
      if strcmp(type1,'TRI3')
%       DKT (FLUIDE, DKIR)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1 1]; % ddl primaux
        nddd = nddp;    % ddl duaux
        nnop = [1 2 3]; % noeuds primaux
        nnod = nnop;    % noeuds duaux
        nnip = nnop;    % fct de base primales
        nnid = nnod;    % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3]; % noeuds pour transformation
        nnit1 = [1 2 3]; % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-2aa')
      end
    elseif strcmp(mode1,'TRID')
%     FLUIDE, TRID
      if (idim ~= 3)
        mode1
        idim
        error('Only 3D for this mode...')
      end
      ddlp = [{'P'}];       % noms des ddl primaux
      ddld = [{'FP'}];      % noms des ddl duaux
      if strcmp(support1,liste_intg{3})
%       FLUIDE, TRID, COMPRESSIBILITE
        comp = [{'P'}];   % noms des composantes primales
        comd = [{'FP'}];  % noms des composantes duales
      elseif strcmp(support1,liste_intg{4})
%       FLUIDE, TRID, PERMEABILITE
        comp = [{'PX'} {'PY'} {'PZ'}];   % noms des composantes primales
        comd = [{'WX'} {'WY'} {'WZ'}];   % noms des composantes duales
      else
        mot1
        mode1
        support1
        error('type of support not implemented 2')
      end
      if strcmp(type1,'SEG2')
%       SEG2 (FLUIDE, TRID)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1]; % ddl primaux
        nddd = nddp;  % ddl duaux
        nnop = [1 2]; % noeuds primaux
        nnod = nnop;  % noeuds duaux
        nnip = nnop;  % fct de base primales
        nnid = nnod;  % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2]; % noeuds pour transformation
        nnit1 = [1 2]; % fct de base pour transformation
      elseif strcmp(type1,'TET4')
%       TET4 (FLUIDE, TRID)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TET4_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TET4_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1 1 1]; % ddl primaux
        nddd = nddp;      % ddl duaux
        nnop = [1 2 3 4]; % noeuds primaux
        nnod = nnop;      % noeuds duaux
        nnip = nnop;      % fct de base primales
        nnid = nnod;      % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4]; % noeuds pour transformation
        nnit1 = [1 2 3 4]; % fct de base pour transformation
      elseif strcmp(type1,'CUB8')
%       CUB8 (FLUIDE, TRID)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = CUB8_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = CUB8_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1 1 1 1 1 1 1]; % ddl primaux
        nddd = nddp;      % ddl duaux
        nnop = [1 2 3 4 5 6 7 8]; % noeuds primaux
        nnod = nnop;      % noeuds duaux
        nnip = nnop;      % fct de base primales
        nnid = nnod;      % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4 5 6 7 8]; % noeuds pour transformation
        nnit1 = [1 2 3 4 5 6 7 8]; % fct de base pour transformation
      elseif strcmp(type1,'PRI6')
%       PRI6 (FLUIDE, TRID)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = PRI6_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = PRI6_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1 1 1 1 1]; % ddl primaux
        nddd = nddp;      % ddl duaux
        nnop = [1 2 3 4 5 6]; % noeuds primaux
        nnod = nnop;      % noeuds duaux
        nnip = nnop;      % fct de base primales
        nnid = nnod;      % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4 5 6]; % noeuds pour transformation
        nnit1 = [1 2 3 4 5 6]; % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-3bis')
      end
    elseif strcmp(mode1,'BARR')
%     FLUIDE, BARR
      ddlp = [{'P'}];       % noms des ddl primaux
      ddld = [{'FP'}];      % noms des ddl duaux
      if strcmp(support1,liste_intg{3})
%       FLUIDE, BARR, COMPRESSIBILITE
        comp = [{'P'}];   % noms des composantes primales
        comd = [{'FP'}];  % noms des composantes duales
      elseif strcmp(support1,liste_intg{4})
%       FLUIDE, BARR, PERMEABILITE
        comp = [{'PX'}];   % noms des composantes primales
        comd = [{'WX'}];   % noms des composantes duales
      else
        mot1
        mode1
        support1
        error('type of support not implemented 2')
      end
      if strcmp(type1,'SEG2')
%       SEG2 (FLUIDE, BARR)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1]; % ddl primaux
        nddd = nddp;  % ddl duaux
        nnop = [1 2]; % noeuds primaux
        nnod = nnop;  % noeuds duaux
        nnip = nnop;  % fct de base primales
        nnid = nnod;  % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2]; % noeuds pour transformation
        nnit1 = [1 2]; % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-3ter')
      end
    else
      mot1
      mode1
      error('mode not implemented 3')
    end

  elseif strcmp(mot1,liste_model{4})

%   THERMIQUE (seuls les noms changent / FLUIDE)
%   """""""""
    if (strcmp(mode1,'COPL') || strcmp(mode1,'DEPL') || strcmp(mode1,'AXIS'))
%     THERMIQUE, (COPL | DEPL | AXIS)
      if (idim ~= 2)
        mode1
        idim
        error('Only 2D for this mode... ter2')
      end
      ddlp = [{'T'}];      % noms des ddl primaux
      ddld = [{'Q'}];      % noms des ddl duaux
      if strcmp(support1,liste_intg{6})
%       THERMIQUE, (COPL | DEPL | AXIS), CAPACITE
        comp = [{'T'}];  % noms des composantes primales
        comd = [{'Q'}];  % noms des composantes duales
      elseif strcmp(support1,liste_intg{7})
%       THERMIQUE, (COPL | DEPL | AXIS), CONDUCTIVITE
        switch mode1
          case {'COPL','DEPL'},
            comp = [{'TX'} {'TY'}];     % noms des composantes primales
            comd = [{'QX'} {'QY'}];     % noms des composantes duales
          case 'AXIS',
            comp = [{'TR'} {'TZ'}];     % noms des composantes primales
            comd = [{'QR'} {'QZ'}];     % noms des composantes duales
          otherwise,
            mode1
            error('Mode not yet implemented THERMIQUE,CONDUCTIVITE... sorry')
        end
      else
        mot1
        mode1
        support1
        error('type of support not implemented')
      end
      if strcmp(type1,list_type_C3M1{2})
%       TRI3 (THERMIQUE, COPL | DEPL)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1 1]; % ddl primaux
        nddd = nddp;    % ddl duaux
        nnop = [1 2 3]; % noeuds primaux
        nnod = nnop;    % noeuds duaux
        nnip = [1 2 3]; % fct de base primales
        nnid = nnip;    % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3];      % noeuds pour transformation
        nnit1 = [1 2 3];      % fct de base pour transformation
      elseif strcmp(type1,list_type_C3M2{2})
%       TRI6 (THERMIQUE, COPL | DEPL)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TRI6_SegmentIntegr2(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TRI6_SegmentIntegr2(support1,mode1);
        end
        nddp = [1 1 1 1 1 1];     % ddl primaux
        nddd = nddp;              % ddl duaux
        nnop = [1 2 3 4 5 6];     % noeuds primaux
        nnod = nnop;              % noeuds duaux
        nnip = [1 2 3 4 5 6];     % fct de base primales
        nnid = nnod;              % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4 5 6];    % noeuds pour transformation
        nnit1 = [1 2 3 4 5 6];    % fct de base pour transformation
      elseif strcmp(type1,list_type_C3M1{3})
%       SEG2 (THERMIQUE, COPL | DEPL)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1]; % ddl primaux
        nddd = nddp;  % ddl duaux
        nnop = [1 2]; % noeuds primaux
        nnod = nnop;  % noeuds duaux
        nnip = nnop;  % fct de base primales
        nnid = nnod;  % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2]; % noeuds pour transformation
        nnit1 = [1 2]; % fct de base pour transformation
      elseif strcmp(type1,'QUA4')
%       QUA4 (THERMIQUE, COPL | DEPL)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = QUA4_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = QUA4_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1 1 1];     % ddl primaux
        nddd = nddp;          % ddl duaux
        nnop = [1 2 3 4];     % noeuds primaux
        nnod = nnop;          % noeuds duaux
        nnip = [1 2 3 4];     % fct de base primales
        nnid = nnod;          % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4];    % noeuds pour transformation
        nnit1 = [1 2 3 4];    % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-3a')
      end
    elseif strcmp(mode1,'TRID')
%     THERMIQUE, TRID
      if (idim ~= 3)
        mode1
        idim
        error('Only 3D for this mode...')
      end
      ddlp = [{'T'}];      % noms des ddl primaux
      ddld = [{'Q'}];      % noms des ddl duaux
      if strcmp(support1,'CAPACITE')
%       THERMIQUE, TRID, CAPACITE
        comp = [{'T'}];  % noms des composantes primales
        comd = [{'Q'}];  % noms des composantes duales
      elseif strcmp(support1,'CONDUCTIVITE')
%       THERMIQUE, TRID, CONDUCTIVITE
        comp = [{'TX'} {'TY'} {'TZ'}];
        comd = [{'QX'} {'QY'} {'QZ'}];
      else
        mot1
        mode1
        support1
        error('type of support not implemented 2a')
      end
      if strcmp(type1,'TET4')
%       TET4 (THERMIQUE, TRID)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TET4_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TET4_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1 1 1]; % ddl primaux
        nddd = nddp;      % ddl duaux
        nnop = [1 2 3 4]; % noeuds primaux
        nnod = nnop;      % noeuds duaux
        nnip = nnop;      % fct de base primales
        nnid = nnod;      % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3 4]; % noeuds pour transformation
        nnit1 = [1 2 3 4]; % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-3ter')
      end
%%DD1
    elseif strcmp(mode1,'BARR')
%     THERMIQUE, BARR
      ddlp = [{'T'}];       % noms des ddl primaux
      ddld = [{'Q'}];       % noms des ddl duaux
      if strcmp(support1,'CAPACITE')
%       THERMIQUE, BARR, CAPACITE
        comp = [{'T'}];   % noms des composantes primales
        comd = [{'Q'}];   % noms des composantes duales
      elseif strcmp(support1,'CONDUCTIVITE')
%       THERMIQUE, BARR, CONDUCTIVITE
        comp = [{'TX'}];   % noms des composantes primales
        comd = [{'QX'}];   % noms des composantes duales
      else
        mot1
        mode1
        support1
        error('type of support not implemented 22')
      end
      if strcmp(type1,'SEG2')
%       SEG2 (THERMIQUE, BARR)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = SEG2_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1]; % ddl primaux
        nddd = nddp;  % ddl duaux
        nnop = [1 2]; % noeuds primaux
        nnod = nnop;  % noeuds duaux
        nnip = nnop;  % fct de base primales
        nnid = nnod;  % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2]; % noeuds pour transformation
        nnit1 = [1 2]; % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-33ter')
      end
%%DD2
    else
      mot1
      mode1
      error('mode not implemented 4')
    end
  elseif strcmp(mot1,liste_model{6})
%   TRANSPORT
    if (strcmp(mode1,liste_mode{2}) || strcmp(mode1,liste_mode{3}))
%     TRANSPORT, (COPL | DEPL)
      ddlp = [{'u'}];      % noms des ddl primaux
      ddld = [{'f'}];      % noms des ddl duaux
      if strcmp(support1,liste_intg{2})
%       TRANSPORT, (COPL | DEPL), MASSE
        comp = [{'u'}];  % noms des composantes primales
        comd = [{'f'}];  % noms des composantes duales
      elseif strcmp(support1,liste_intg{8})
%       TRANSPORT, (COPL | DEPL), ADVECTION
        switch idim
          case 2,
            comp = [{'uX'} {'uY'}];     % noms des composantes primales
            comd = [{'fX'} {'fY'}];     % noms des composantes duales
%          case 3,
%            comp = [{'TX'} {'TY'} {'TZ'}];
%            comd = [{'QX'} {'QY'} {'QZ'}];
          otherwise,
            error('bad idim')
        end
      elseif strcmp(support1,liste_intg{9})
%       TRANSPORT, (COPL | DEPL), FLUX
        switch idim
          case 2,
            comp = [{'U'}];     % noms des composantes primales
            comd = [{'F'}];     % noms des composantes duales
          otherwise,
            error('bad idim')
        end
      else
        mot1
        mode1
        support1
        error('type of support not implemented-7')
      end
      if strcmp(type1,list_type_C3M1{2})
%       TRI3 (TRANSPORT, COPL | DEPL)
        if exist('intg2')
          [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1, ...
                                                       intg2{zo1});
        else
          [PHI,DPHI,COOR,WEIGHT] = TRI3_SegmentIntegr3(support1,mode1);
        end
        nddp = [1 1 1]; % ddl primaux
        nddd = nddp;    % ddl duaux
        nnop = [1 2 3]; % noeuds primaux
        nnod = nnop;    % noeuds duaux
        nnip = [1 2 3]; % fct de base primales
        nnid = nnip;    % fct de base duales
        ncop = [1:length(comp)];  % numeros des composantes primales
        ncod = [1:length(comd)];  % numeros des composantes duales
        nnot1 = [1 2 3];      % noeuds pour transformation
        nnit1 = [1 2 3];      % fct de base pour transformation
      else
        mot1
        mode1
        type1
        error('type of element not implemented-7')
      end
    else
      mot1
      mode1
      error('mode not implemented 7')
    end
  else
    mot1
    error('model not implemented')
  end
%
  modl1{zo1} = struct('DDLP',{ddlp},'DDLD',{ddld}, ...
                      'COMP',{comp},'COMD',{comd}, ...
                      'NNOP',nnop,'NNOD',nnod, ...
		      'NDDP',nddp,'NDDD',nddd, ...
		      'NNIP',nnip,'NNID',nnid, ...
                      'NCOP',ncop,'NCOD',ncod, ...
		      'NNOT',nnot1, ...
		      'NNIT',nnit1);
% For elements of structures (beam, plate...),
% local basis dof
  if exist('ddlploc')
    modl1{zo1}.DDLPLOC = ddlploc;
    modl1{zo1}.NDDPLOC = nddploc;
    modl1{zo1}.NNOPLOC = nnoploc;
  end
  if exist('ddldloc')
    modl1{zo1}.DDLDLOC = ddldloc;
    modl1{zo1}.NDDDLOC = ndddloc;
    modl1{zo1}.NNODLOC = nnodloc;
  end
  intg1{zo1} = struct('PHI',PHI,'DPHI',DPHI, ...
                      'COOR',COOR,'WEIGHT',WEIGHT);
% For elements of structures (beam, plate...),
% shape function derivates at nodes
  if exist('DPHIN')
    intg1{zo1}.DPHIN = DPHIN;
  end
  if exist('DDPHI')
    intg1{zo1}.DDPHI = DDPHI;
  end
end
