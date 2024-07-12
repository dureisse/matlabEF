function [dshp1] = EF_Dshape(type1,X);
% First derivative of finite element shape functions 
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 23 / 01 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 12 / 04 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 21 / 07 / 2003
%   Ajout SEG2
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 15 / 07 / 2004
%   Ajout DKQ, DKT
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 13 / 12 / 2004
%   Ajout CUB8
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 15 / 05 / 2005
%   Ajout RAC2
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 16 / 09 / 2005
%   Ajout RAC3, QUA8
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 02 / 07 / 2006
%   Ajout TET4, PRI6
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 15 / 12 / 2006
%   Ajout POI1
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 03 / 03 / 2007
%   Ajout SEG3, CU20, PR15
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 24 / 11 / 2007
%   Ajout du TET5b, TE10
% DUREISSEIX David  LaMCoS                          le 08 / 10 / 2010
%   Ajout DSQ, DST
% DUREISSEIX David  LaMCoS                          le 15 / 06 / 2015
%   correction de bug sur QUA8
% 
% Retourne les valeurs des derivees des fonctions de forme d'un element
% geometrique de type type1 au point dont les coordonnees de reference
% sont X
%   
% Entrees
%   type1               : type de l'element geometrique
%   X(idimr)            : coordonnees de reference du point
% Sorties
%   dshp1(idimr,nbnn)   : valeurs des derivees en ce point

switch type1
  case 'POI1',
    dshp1 = 0.; % par convention

  case 'SEG2',
    dshp1 = 0.5 * [-1. 1.];

  case 'SEG3',
    dshp1 = [(X(1)-0.5) (X(1)+0.5) -2.*X(1)];

  case 'TRI3',
    dshp1 = [-1. 1. 0.
             -1. 0. 1.];

  case 'QUA4',
    dshp1 = 0.25 * [ (X(2)-1.) (1.-X(2)) (1.+X(2)) -(1.+X(2))
                     (X(1)-1.) -(1.+X(1)) (1.+X(1)) (1.-X(1))];

  case 'TRI6',
    dshp1 = [4.*X(2)+4.*X(1)-3.     4.*X(2)+4.*X(1)-3.
             4.*X(1)-1.             0.
             0.                     4.*X(2)-1.
             -4.*(X(2)+2.*X(1)-1.)  -4.*X(1)
             4.*X(2)                4.*X(1)
             -4.*X(2)               -4.*(2.*X(2)+X(1)-1.)]';

  case {'DKQ','DSQ'}
    dshp1 = 0.25 * [ (X(2)-1.) (1.-X(2)) (1.+X(2)) -(1.+X(2))
                     (X(1)-1.) -(1.+X(1)) (1.+X(1)) (1.-X(1))];
    dshp2 = [-X(1)*(1.-X(2)) 0.5*(1.-X(2)^2) -X(1)*(1.+X(2)) -0.5*(1.-X(2)^2)
             -0.5*(1.-X(1)^2) -(1.+X(1))*X(2) 0.5*(1.-X(1)^2) -(1.-X(1))*X(2)];
    dshp1 = [dshp1 dshp2];
    clear dshp2;

  case {'DKT','DST'}
    dshp1 = [-1. 1. 0.
             -1. 0. 1.];
    dshp2 = 4. * [(1.-2.*X(1)-X(2)) X(2) -X(2)
                  -X(1)             X(1) (1.-X(1)-2.*X(2))];
    dshp1 = [dshp1 dshp2];
    clear dshp2;

  case 'CUB8',
%   dshp1(idimr,nbnn)   : valeurs des derivees en ce point
    QSIM= 1-X(1); QSIP= 1+X(1);
    ETAM= 1-X(2); ETAP= 1+X(2);
    DZEM= 1-X(3); DZEP= 1+X(3);
    dshp1 = [-ETAM*DZEM/8 -QSIM*DZEM/8 -QSIM*ETAM/8
              ETAM*DZEM/8 -QSIP*DZEM/8 -QSIP*ETAM/8
              ETAP*DZEM/8  QSIP*DZEM/8 -QSIP*ETAP/8
             -ETAP*DZEM/8  QSIM*DZEM/8 -QSIM*ETAP/8
             -ETAM*DZEP/8 -QSIM*DZEP/8  QSIM*ETAM/8
              ETAM*DZEP/8 -QSIP*DZEP/8  QSIP*ETAM/8
              ETAP*DZEP/8  QSIP*DZEP/8  QSIP*ETAP/8
             -ETAP*DZEP/8  QSIM*DZEP/8  QSIM*ETAP/8]';

  case 'CU20',
%   D'apres cast3m shape3.eso
      UNQUA = 0.25; UNDEMI = 0.5; UN = 1.; DEUX = 2.; CINQ = 5.; HUIT = 8.;
      QSI = X(1);
      ETA = X(2);
      DZE = X(3);
      QSIM=(UN-QSI);
      QSIP=(UN+QSI);
      ETAM=(UN-ETA);
      ETAP=(UN+ETA);
      DZEM=(UN-DZE);
      DZEP=(UN+DZE);
      SHP = zeros(4,20);
%     NOEUDS AUX SOMMETS 1 2 3 4 5 6 7 8
      SHP(2,1 )=    -ETAM*DZEM*(DEUX*QSIM+ETAM+DZEM-CINQ)/HUIT;
      SHP(2,2 )=     ETAM*DZEM*(DEUX*QSIP+ETAM+DZEM-CINQ)/HUIT;
      SHP(2,3 )=     ETAP*DZEM*(DEUX*QSIP+ETAP+DZEM-CINQ)/HUIT;
      SHP(2,4 )=    -ETAP*DZEM*(DEUX*QSIM+ETAP+DZEM-CINQ)/HUIT;
      SHP(2,5 )=    -ETAM*DZEP*(DEUX*QSIM+ETAM+DZEP-CINQ)/HUIT;
      SHP(2,6 )=     ETAM*DZEP*(DEUX*QSIP+ETAM+DZEP-CINQ)/HUIT;
      SHP(2,7 )=     ETAP*DZEP*(DEUX*QSIP+ETAP+DZEP-CINQ)/HUIT;
      SHP(2,8 )=    -ETAP*DZEP*(DEUX*QSIM+ETAP+DZEP-CINQ)/HUIT;
      SHP(3,1 )=    -QSIM*DZEM*(QSIM+DEUX*ETAM+DZEM-CINQ)/HUIT;
      SHP(3,2 )=    -QSIP*DZEM*(QSIP+DEUX*ETAM+DZEM-CINQ)/HUIT;
      SHP(3,3 )=     QSIP*DZEM*(QSIP+DEUX*ETAP+DZEM-CINQ)/HUIT;
      SHP(3,4 )=     QSIM*DZEM*(QSIM+DEUX*ETAP+DZEM-CINQ)/HUIT;
      SHP(3,5 )=    -QSIM*DZEP*(QSIM+DEUX*ETAM+DZEP-CINQ)/HUIT;
      SHP(3,6 )=    -QSIP*DZEP*(QSIP+DEUX*ETAM+DZEP-CINQ)/HUIT;
      SHP(3,7 )=     QSIP*DZEP*(QSIP+DEUX*ETAP+DZEP-CINQ)/HUIT;
      SHP(3,8 )=     QSIM*DZEP*(QSIM+DEUX*ETAP+DZEP-CINQ)/HUIT;
      SHP(4,1 )=    -QSIM*ETAM*(QSIM+ETAM+DEUX*DZEM-CINQ)/HUIT;
      SHP(4,2 )=    -QSIP*ETAM*(QSIP+ETAM+DEUX*DZEM-CINQ)/HUIT;
      SHP(4,3 )=    -QSIP*ETAP*(QSIP+ETAP+DEUX*DZEM-CINQ)/HUIT;
      SHP(4,4 )=    -QSIM*ETAP*(QSIM+ETAP+DEUX*DZEM-CINQ)/HUIT;
      SHP(4,5 )=     QSIM*ETAM*(QSIM+ETAM+DEUX*DZEP-CINQ)/HUIT;
      SHP(4,6 )=     QSIP*ETAM*(QSIP+ETAM+DEUX*DZEP-CINQ)/HUIT;
      SHP(4,7 )=     QSIP*ETAP*(QSIP+ETAP+DEUX*DZEP-CINQ)/HUIT;
      SHP(4,8 )=     QSIM*ETAP*(QSIM+ETAP+DEUX*DZEP-CINQ)/HUIT;
%     NOEUDS SUR LES COTES PARALLELES A L AXE QSI 9 11 13 15
      SHP(2,9 )=-UNDEMI*QSI*ETAM*DZEM;
      SHP(2,11)=-UNDEMI*QSI*ETAP*DZEM;
      SHP(2,13)=-UNDEMI*QSI*ETAM*DZEP;
      SHP(2,15)=-UNDEMI*QSI*ETAP*DZEP;
      SHP(3,9 )=-UNQUA*(UN-QSI*QSI)*DZEM;
      SHP(3,11)= UNQUA*(UN-QSI*QSI)*DZEM;
      SHP(3,13)=-UNQUA*(UN-QSI*QSI)*DZEP;
      SHP(3,15)= UNQUA*(UN-QSI*QSI)*DZEP;
      SHP(4,9 )=-UNQUA*(UN-QSI*QSI)*ETAM;
      SHP(4,11)=-UNQUA*(UN-QSI*QSI)*ETAP;
      SHP(4,13)= UNQUA*(UN-QSI*QSI)*ETAM;
      SHP(4,15)= UNQUA*(UN-QSI*QSI)*ETAP;
%     NOEUDS SUR LES COTES PARALELLES A L AXE ETA 10 12 14 16
      SHP(2,10)= UNQUA*(UN-ETA*ETA)*DZEM;
      SHP(2,12)=-UNQUA*(UN-ETA*ETA)*DZEM;
      SHP(2,14)= UNQUA*(UN-ETA*ETA)*DZEP;
      SHP(2,16)=-UNQUA*(UN-ETA*ETA)*DZEP;
      SHP(3,10)=-UNDEMI*QSIP*ETA*DZEM;
      SHP(3,12)=-UNDEMI*QSIM*ETA*DZEM;
      SHP(3,14)=-UNDEMI*QSIP*ETA*DZEP;
      SHP(3,16)=-UNDEMI*QSIM*ETA*DZEP;
      SHP(4,10)=-UNQUA*QSIP*(UN-ETA*ETA);
      SHP(4,12)=-UNQUA*QSIM*(UN-ETA*ETA);
      SHP(4,14)= UNQUA*QSIP*(UN-ETA*ETA);
      SHP(4,16)= UNQUA*QSIM*(UN-ETA*ETA);
%     NOEUDS SUR LES COTES PARALELLES A L AXE DZE 17 18 19 20
      SHP(2,17)=-UNQUA*ETAM*(UN-DZE*DZE);
      SHP(2,18)= UNQUA*ETAM*(UN-DZE*DZE);
      SHP(2,19)= UNQUA*ETAP*(UN-DZE*DZE);
      SHP(2,20)=-UNQUA*ETAP*(UN-DZE*DZE);
      SHP(3,17)=-UNQUA*QSIM*(UN-DZE*DZE);
      SHP(3,18)=-UNQUA*QSIP*(UN-DZE*DZE);
      SHP(3,19)= UNQUA*QSIP*(UN-DZE*DZE);
      SHP(3,20)= UNQUA*QSIM*(UN-DZE*DZE);
      SHP(4,17)=-UNDEMI*QSIM*ETAM*DZE;
      SHP(4,18)=-UNDEMI*QSIP*ETAM*DZE;
      SHP(4,19)=-UNDEMI*QSIP*ETAP*DZE;
      SHP(4,20)=-UNDEMI*QSIM*ETAP*DZE;
    dshp1 = SHP(2:4,:);
    clear SHP;

  case 'RAC2',
    dshp1 = [-1. 1.];

  case 'RAC3',
    disp('WarningEFDshape : RAC3 non teste (node order?)')
    dshp1 = [(X(1)-0.5) -2.*X(1) (X(1)+0.5)];

  case 'QUA8',
    QSIM= 1-X(1); QSIP= 1+X(1);
    ETAM= 1-X(2); ETAP= 1+X(2);
%     dshp1 = [ETAM*(2*X(1)+X(2))/4    QSIM*(X(1)+2*X(2))/4
%              -X(1)*ETAM              QSIM*ETAM/-2 !!! BUG cf ci-apres
% 	         ETAM*(2*X(1)-X(2))/4    QSIP*(X(1)-2*X(2))/-4
%              ETAP*ETAM/2             -X(2)*QSIP
%              ETAP*(2*X(1)+X(2))/4    QSIP*(X(1)+2*X(2))/4
% 	         -X(1)*ETAP              QSIP*QSIM/2
%              ETAP*(2*X(1)-X(2))/4    QSIM*(X(1)-2*X(2))/-4
% 	         ETAP*ETAM/-2            -X(2)*QSIM]';
    dshp1 = [ETAM*(2*X(1)+X(2))/4    QSIM*(X(1)+2*X(2))/4
	         ETAM*(2*X(1)-X(2))/4    QSIP*(X(1)-2*X(2))/-4
             ETAP*(2*X(1)+X(2))/4    QSIP*(X(1)+2*X(2))/4
             ETAP*(2*X(1)-X(2))/4    QSIM*(X(1)-2*X(2))/-4
             -X(1)*ETAM              QSIM*QSIP/-2
             ETAP*ETAM/2             -X(2)*QSIP
	         -X(1)*ETAP              QSIP*QSIM/2
	         ETAP*ETAM/-2            -X(2)*QSIM]';

  case 'TET4',
%   dshp1(idimr,nbnn)   : valeurs des derivees en ce point
    dshp1          = [-1. 1. 0. 0.
                      -1. 0. 1. 0.
		      -1. 0. 0. 1.];
  case 'TET5b',
%   dshp1(idimr,nbnn)   : valeurs des derivees en ce point
    disp('Warning: TET5n on test only (not fixed yet)')
% TET5b is a TET4 element plus a bubble-like 5th node at the centroid
% The "bubble" is a pyramid-like shape function (linear on each
% sub-triangle), non differentiable on the radial lines (!)
% Original shape functions of the TET4 element are modified to
% preserve interpolation property, see Hughes (The FEM, linear
% static and dynamic FEA).
%   Shape functions of TET4 = barycentric coordinates
    shp0 = [(1.-X(1)-X(2)-X(3)) X(1) X(2) X(3)];
%   Derivative of shape functions of TET4
    dshp0          = [-1. 1. 0. 0.
                      -1. 0. 1. 0.
		      -1. 0. 0. 1.];
%   Which node is the less close? This gives the 5th shape function
    [phii,i] = min(shp0);
    dphi5 = 4.*dshp0(:,i);
%   Corrected TET4 derivatives of shape functions
%   and derivatives of 5th shape function
    dshp1 = [dshp0 - 0.25*kron([1. 1. 1. 1.],dphi5) , dphi5];

  case 'TE10',
%     D'apres cast3m shape3.eso
      UNQUA = 0.25; UNDEMI = 0.5; UN = 1.; DEUX = 2.; CINQ = 5.; HUIT = 8.;
      QUATRE = 4.; XZER = 0.;
      QSI = X(1);
      ETA = X(2);
      DZE = X(3);
      QSIM=(UN-QSI);
      QSIP=(UN+QSI);
      ETAM=(UN-ETA);
      ETAP=(UN+ETA);
      DZEM=(UN-DZE);
      DZEP=(UN+DZE);
      SHP = zeros(4,10);
      AUX=UN-QSI-ETA-DZE;
%     NOEUDS AUX SOMMETS 1 2 3 4 (dans cast3m : 1 3 5 10)
      SHP(2,1)= UN-QUATRE*AUX;
      SHP(3,1)= UN-QUATRE*AUX;
      SHP(4,1)= UN-QUATRE*AUX;
      SHP(2,2)= QUATRE*QSI-UN;
      SHP(3,2)= XZER;
      SHP(4,2)= XZER;
      SHP(2,3)= XZER;
      SHP(3,3)= QUATRE*ETA-UN;
      SHP(4,3)= XZER;
      SHP(2,4)= XZER;
      SHP(3,4)= XZER;
      SHP(4,4)= QUATRE*DZE-UN;
%     NOEUD MILIEU ARETE 1 2 (dans cast3m 2)
      SHP(2,5)= QUATRE*(AUX-QSI);
      SHP(3,5)= -QUATRE*QSI;
      SHP(4,5)= -QUATRE*QSI;
%     NOEUD MILIEU ARETE 2 3 (dans cast3m 4)
      SHP(2,6)= QUATRE*ETA;
      SHP(3,6)= QUATRE*QSI;
      SHP(4,6)= XZER;
%     NOEUD MILIEU ARETE 3 1 (dans cast3m 6)
      SHP(2,7)= -QUATRE*ETA;
      SHP(3,7)= QUATRE*(AUX-ETA);
      SHP(4,7)= -QUATRE*ETA;
%     NOEUD MILIEU ARETE 4 1 (dans cast3m 7)
      SHP(2,8)= -QUATRE*DZE;
      SHP(3,8)= -QUATRE*DZE;
      SHP(4,8)= QUATRE*(AUX-DZE);
%     NOEUD MILIEU ARETE 4 3 (dans cast3m 9)
      SHP(2,9)= XZER;
      SHP(3,9)= QUATRE*DZE;
      SHP(4,9)= QUATRE*ETA;
%     NOEUD MILIEU ARETE 4 2 (dans cast3m 8)
      SHP(2,10)= QUATRE*DZE;
      SHP(3,10)= XZER;
      SHP(4,10)= QUATRE*QSI;
    dshp1 = SHP(2:4,:);
    clear SHP;

  case 'PRI6',
%   dshp1(idimr,nbnn)   : valeurs des derivees en ce point
    dshp1(1:3,1:6) = [-0.5*(1.-X(3)) -0.5*(1.-X(3)) -0.5*(1.-X(1)-X(2))
                       0.5*(1.-X(3))  0.            -0.5*X(1)
		       0.             0.5*(1.-X(3)) -0.5*X(2)
		      -0.5*(1.+X(3)) -0.5*(1.+X(3))  0.5*(1.-X(1)-X(2))
		       0.5*(1.+X(3))  0.             0.5*X(1)
		       0.             0.5*(1.+X(3))  0.5*X(2)]';

  case 'PR15',
%   dshp1(idimr,nbnn)   : valeurs des derivees en ce point
%   D'apres cast3m shape3.eso
      XZER = 0.; UNQUA = 0.25; UNDEMI = 0.5; UN = 1.; DEUX = 2.; CINQ = 5.; HUIT = 8.;
      QSI = X(1);
      ETA = X(2);
      DZE = X(3);
      QSIM=(UN-QSI);
      QSIP=(UN+QSI);
      ETAM=(UN-ETA);
      ETAP=(UN+ETA);
      DZEM=(UN-DZE);
      DZEP=(UN+DZE);
      SHP = zeros(4,15);
      AUX=UN-QSI-ETA;
      PAUX =QSI-UNDEMI*DZE-UN;
      PAUX1=ETA-UNDEMI*DZE-UN;
      PAUX2=QSI+ETA+UNDEMI*DZE;
      PAUX3=QSI+UNDEMI*DZE-UN;
      PAUX4=ETA+UNDEMI*DZE-UN;
      PAUX5=QSI+ETA-UNDEMI*DZE;
      SHP(2,1 )=(PAUX2-AUX)*DZEM;
      SHP(2,7 )=DEUX*(AUX-QSI)*DZEM;
      SHP(2,2 )=(PAUX+QSI)*DZEM;
      SHP(3,1 )=(PAUX2-AUX)*DZEM;
      SHP(3,7 )=-DEUX*QSI*DZEM;
      SHP(3,2 )=XZER;
      SHP(4,1 )=AUX*(PAUX2-UNDEMI*DZEM);
      SHP(4,7 )=-DEUX*QSI*AUX;
      SHP(4,2 )=-QSI*(UNDEMI*DZEM+PAUX);
      SHP(2,8 )=DEUX*ETA*DZEM;
      SHP(2,3 )=XZER;
      SHP(2,9 )=-DEUX*ETA*DZEM;
      SHP(3,8 )=DEUX*QSI*DZEM;
      SHP(3,3 )=(ETA+PAUX1)*DZEM;
      SHP(3,9 )=DEUX*(AUX-ETA)*DZEM;
      SHP(4,8 )=-DEUX*QSI*ETA;
      SHP(4,3 )=-ETA*(PAUX1+UNDEMI*DZEM);
      SHP(4,9 )=-DEUX*ETA*AUX;
      SHP(2,13)=-DZEM*DZEP;
      SHP(2,14)=DZEM*DZEP;
      SHP(2,15)=XZER;
      SHP(3,13)=-DZEM*DZEP;
      SHP(3,14)=XZER;
      SHP(3,15)=DZEM*DZEP;
      SHP(4,13)=-DEUX*DZE*AUX;
      SHP(4,14)=-DEUX*DZE*QSI;
      SHP(4,15)=-DEUX*DZE*ETA;
      SHP(2,4 )=(PAUX5-AUX)*DZEP;
      SHP(2,10)=DEUX*(AUX-QSI)*DZEP;
      SHP(2,5 )=(PAUX3+QSI)*DZEP;
      SHP(3,4 )=(PAUX5-AUX)*DZEP;
      SHP(3,10)=-DEUX*QSI*DZEP;
      SHP(3,5 )=XZER;
      SHP(4,4 )=AUX*(UNDEMI*DZEP-PAUX5);
      SHP(4,10)=DEUX*QSI*AUX;
      SHP(4,5 )=QSI*(PAUX3+UNDEMI*DZEP);
      SHP(2,11)=DEUX*ETA*DZEP;
      SHP(2,6 )=XZER;
      SHP(2,12)=-DEUX*ETA*DZEP;
      SHP(3,11)=DEUX*QSI*DZEP;
      SHP(3,6 )=(ETA+PAUX4)*DZEP;
      SHP(3,12)=DEUX*(AUX-ETA)*DZEP;
      SHP(4,11)=DEUX*QSI*ETA;
      SHP(4,6 )=ETA*(PAUX4+UNDEMI*DZEP);
      SHP(4,12)=DEUX*ETA*AUX;
    dshp1 = SHP(2:4,:);
    clear SHP;

  otherwise,
    type1
    error('Type of element not implemented yet')
end

