function [shp1] = EF_Shape(type1,X);
% Values of finite element shape functions
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 23 / 01 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 12 / 04 / 2003
%   Ajout du TRI6
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 21 / 07 / 2003
%   Ajout du SEG2
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 12 / 07 / 2004
%   Ajout du DKQ, DKT
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 13 / 12 / 2004
%   Ajout du CUB8
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 15 / 05 / 2005
%   Ajout du RAC2
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 16 / 09 / 2005
%   Ajout du RAC3, QUA8
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 02 / 07 / 2006
%   Ajout du TET4, PRI6
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 15 / 12 / 2006
%   Ajout du POI1
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 03 / 03 / 2007
%   Ajout du SEG3, CU20
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 24 / 11 / 2007
%   Ajout du TET5b, TE10
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 09 / 08 / 2008
%   Ajout du TET9b
% DUREISSEIX David  LaMCoS                          le 08 / 10 / 2010
%   Ajout du DSQ,DST
% 
% Retourne les valeurs des fonctions de forme d'un element geometrique
% de type type1 au point dont les coordonnees de reference sont X
%   
% Entrees
%   type1               : type de l'element geometrique
%   X(idimr)            : coordonnees de reference du point
% Sorties
%   shp1(nbnn)          : valeurs des fonctions de forme en ce point

switch type1
  case 'POI1',
    shp1 = 1.; % par convention

  case 'SEG2',
    shp1 = 0.5*[1.-X(1)   1.+X(1)];

  case 'SEG3',
    shp1 = [0.5*X(1)*(X(1)-1.) 0.5*X(1)*(X(1)+1.) (1.-X(1)*X(1))];

  case 'TRI3',
    shp1 = [(1.-X(1)-X(2))   X(1)   X(2)];

  case 'QUA4',
    shp1 = 0.25 * [(1.-X(1))*(1.-X(2)) (1.+X(1))*(1.-X(2)) ...
                   (1.+X(1))*(1.+X(2)) (1.-X(1))*(1.+X(2))];

  case 'TRI6',
    shp1 = [(2*X(2)+2.*X(1)-1.)*(X(2)+X(1)-1.)
            X(1)*(2.*X(1)-1.)
            X(2)*(2.*X(2)-1.)
            -4.*X(1)*(X(2)+X(1)-1.)
            4.*X(1)*X(2)
            -4.*X(2)*(X(2)+X(1)-1.)]';

  case {'DKQ','DSQ'}
    shp1 = 0.25 * [(1.-X(1))*(1.-X(2)) (1.+X(1))*(1.-X(2)) ...
                   (1.+X(1))*(1.+X(2)) (1.-X(1))*(1.+X(2))];
    shp2 = 0.5 * [(1.-X(1)^2)*(1.-X(2)) (1.+X(1))*(1.-X(2)^2) ...
                  (1.-X(1)^2)*(1.+X(2)) (1.-X(1))*(1.-X(2)^2)];
    shp1 = [shp1 shp2];
    clear shp2;

  case {'DKT','DST'}
    shp1 = [(1.-X(1)-X(2)) X(1) X(2)];
    shp2 = 4. * [X(1)*(1.-X(1)-X(2)) X(1)*X(2) X(2)*(1.-X(1)-X(2))];
    shp1 = [shp1 shp2];
    clear shp2;

  case 'CUB8',
    QSIM= 1-X(1); QSIP= 1+X(1);
    ETAM= 1-X(2); ETAP= 1+X(2);
    DZEM= 1-X(3); DZEP= 1+X(3);
    shp1 = [QSIM*ETAM*DZEM/8
            QSIP*ETAM*DZEM/8
            QSIP*ETAP*DZEM/8
            QSIM*ETAP*DZEM/8
            QSIM*ETAM*DZEP/8
            QSIP*ETAM*DZEP/8
            QSIP*ETAP*DZEP/8
            QSIM*ETAP*DZEP/8]';

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
      SHP = zeros(1,20);
%     NOEUDS AUX SOMMETS 1 2 3 4 5 6 7 8
      SHP(1,1 )=QSIM*ETAM*DZEM*(QSIM+ETAM+DZEM-CINQ)/HUIT;
      SHP(1,2 )=QSIP*ETAM*DZEM*(QSIP+ETAM+DZEM-CINQ)/HUIT;
      SHP(1,3 )=QSIP*ETAP*DZEM*(QSIP+ETAP+DZEM-CINQ)/HUIT;
      SHP(1,4 )=QSIM*ETAP*DZEM*(QSIM+ETAP+DZEM-CINQ)/HUIT;
      SHP(1,5 )=QSIM*ETAM*DZEP*(QSIM+ETAM+DZEP-CINQ)/HUIT;
      SHP(1,6 )=QSIP*ETAM*DZEP*(QSIP+ETAM+DZEP-CINQ)/HUIT;
      SHP(1,7 )=QSIP*ETAP*DZEP*(QSIP+ETAP+DZEP-CINQ)/HUIT;
      SHP(1,8 )=QSIM*ETAP*DZEP*(QSIM+ETAP+DZEP-CINQ)/HUIT;
%     NOEUDS SUR LES COTES PARALLELES A L AXE QSI 9 11 13 15
      SHP(1,9 )=UNQUA*(UN-QSI*QSI)*ETAM*DZEM;
      SHP(1,11)=UNQUA*(UN-QSI*QSI)*ETAP*DZEM;
      SHP(1,13)=UNQUA*(UN-QSI*QSI)*ETAM*DZEP;
      SHP(1,15)=UNQUA*(UN-QSI*QSI)*ETAP*DZEP;
%     NOEUDS SUR LES COTES PARALELLES A L AXE ETA 10 12 14 16
      SHP(1,10)=UNQUA*QSIP*(UN-ETA*ETA)*DZEM;
      SHP(1,12)=UNQUA*QSIM*(UN-ETA*ETA)*DZEM;
      SHP(1,14)=UNQUA*QSIP*(UN-ETA*ETA)*DZEP;
      SHP(1,16)=UNQUA*QSIM*(UN-ETA*ETA)*DZEP;
%     NOEUDS SUR LES COTES PARALELLES A L AXE DZE 17 18 19 20
      SHP(1,17)= UNQUA*QSIM*ETAM*(UN-DZE*DZE);
      SHP(1,18)= UNQUA*QSIP*ETAM*(UN-DZE*DZE);
      SHP(1,19)= UNQUA*QSIP*ETAP*(UN-DZE*DZE);
      SHP(1,20)= UNQUA*QSIM*ETAP*(UN-DZE*DZE);
    shp1 = SHP(1,:);
    clear SHP;

  case 'RAC2',
    shp1 = [1.-X(1)   X(1)];

  case 'RAC3',
    disp('Warning: RAC3 non teste (node order?)')
    shp1 = [-0.5*(1.-X(1))*X(1) (1.-X(1))*(1.+X(1)) 0.5*(1.+X(1))*X(1)];

  case 'QUA8',
    QSIM= 1-X(1); QSIP= 1+X(1);
    ETAM= 1-X(2); ETAP= 1+X(2);
%     shp1 = [QSIM*ETAM*(1+X(1)+X(2))/-4
%             QSIP*QSIM*ETAM/2
% 	        QSIP*ETAM*(1-X(1)+X(2))/-4
% 	        QSIP*ETAP*ETAM/2
% 	        QSIP*ETAP*(-1+X(1)+X(2))/4
% 	        QSIP*QSIM*ETAP/2
% 	        QSIM*ETAP*(1+X(1)-X(2))/-4
% 	        QSIM*ETAM*ETAP/2]';
    shp1 = [QSIM*ETAM*(1+X(1)+X(2))/-4
	        QSIP*ETAM*(1-X(1)+X(2))/-4
	        QSIP*ETAP*(-1+X(1)+X(2))/4
	        QSIM*ETAP*(1+X(1)-X(2))/-4
            QSIP*QSIM*ETAM/2
	        QSIP*ETAP*ETAM/2
	        QSIP*QSIM*ETAP/2
	        QSIM*ETAM*ETAP/2]';

  case 'TET4',
    shp1 = [(1.-X(1)-X(2)-X(3)) X(1) X(2) X(3)];

  case 'TET5b',
    disp('Warning: TET5b on test only (not fixed yet)')
% TET5b is a TET4 element plus a bubble-like 5th node at the centroid
% The "bubble" is a pyramid-like shape function (linear on each
% sub-triangle), non differentiable on the radial lines (!)
% Original shape functions of the TET4 element are modified to 
% preserve interpolation property, see Hughes (The FEM, linear
% static and dynamic FEA).
%   Shape functions of TET4 = barycentric coordinates
    shp0 = [(1.-X(1)-X(2)-X(3)) X(1) X(2) X(3)];
%   Which node is the less close? This gives the 5th shape function
    phi5 = 4.*min(shp0);
%   Corrected TET4 shape functions to recover interpolation property
%   and 5th shape function
    shp1 = [shp0 - 0.25*phi5 , phi5];

  case 'TET9b',
    disp('Warning: TET9b on test only (not fixed yet)')
    error('TO BE DONE')
% TET9b is a TET4 element plus a bubble-like defined at 4 additional
% internal nodes. It is build recursively with TET5b applied to the
% TET4 element, and each sub TET4 is again transformed into TET5b.
%   Shape functions of TET4 = barycentric coordinates
    shp0 = [(1.-X(1)-X(2)-X(3)) X(1) X(2) X(3)];
%   Which node is the less close? This gives the 5th shape function
    phi5 = 4.*min(shp0);
%   Corrected TET4 shape functions to recover interpolation property
%   and 5th shape function
    shp5 = [shp0 - 0.25*phi5 , phi5];
%   Again...

  case 'TE10',
%     D'apres cast3m shape3.eso
      UNQUA = 0.25; UNDEMI = 0.5; UN = 1.; DEUX = 2.; CINQ = 5.; HUIT = 8.;
      QUATRE = 4.;
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
%     NOEUDS AUX SOMMETS 1 2 3 4
%     SHP(1,1)=AUX*(DEUX*AUX-UN)
      SHP(1,1)=AUX*(DEUX*AUX-UN);
%     SHP(1,3)=QSI*(DEUX*QSI-UN)
      SHP(1,2)=QSI*(DEUX*QSI-UN);
%     SHP(1,5)=ETA*(DEUX*ETA-UN)
      SHP(1,3)=ETA*(DEUX*ETA-UN);
%     SHP(1,10)=DZE*(DEUX*DZE-UN)
      SHP(1,4)=DZE*(DEUX*DZE-UN);
%     NOEUD MILIEU ARETE 1 2
%     SHP(1,2)=QUATRE*QSI*AUX
      SHP(1,5)=QUATRE*QSI*AUX;
%     NOEUD MILIEU ARETE 2 3
%     SHP(1,4)=QUATRE*QSI*ETA
      SHP(1,6)=QUATRE*QSI*ETA;
%     NOEUD MILIEU ARETE 3 1
%     SHP(1,6)=QUATRE*ETA*AUX
      SHP(1,7)=QUATRE*ETA*AUX;
%     NOEUD MILIEU ARETE 4 1
%     SHP(1,7)=QUATRE*DZE*AUX
      SHP(1,8)=QUATRE*DZE*AUX;
%     NOEUD MILIEU ARETE 4 3
%     SHP(1,9)=QUATRE*ETA*DZE
      SHP(1,9)=QUATRE*ETA*DZE;
%     NOEUD MILIEU ARETE 4 2
%     SHP(1,8)=QUATRE*QSI*DZE
      SHP(1,10)=QUATRE*QSI*DZE;
    shp1 = SHP(1,:);
    clear SHP;

  case 'PRI6',
    shp1 = [0.5*(1.-X(1)-X(2))*(1.-X(3)) ...
            0.5*X(1)*(1.-X(3)) ...
            0.5*X(2)*(1.-X(3)) ...
	    0.5*(1.-X(1)-X(2))*(1.+X(3)) ...
	    0.5*X(1)*(1.+X(3)) ...
	    0.5*X(2)*(1.+X(3))];

  case 'PR15',
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
      SHP = zeros(4,15);
      AUX=UN-QSI-ETA;
      PAUX =QSI-UNDEMI*DZE-UN;
      PAUX1=ETA-UNDEMI*DZE-UN;
      PAUX2=QSI+ETA+UNDEMI*DZE;
      PAUX3=QSI+UNDEMI*DZE-UN;
      PAUX4=ETA+UNDEMI*DZE-UN;
      PAUX5=QSI+ETA-UNDEMI*DZE;
%     ORDRE DES NOEUDS DE cast3M : 1 2 3 4 5 6 7  8  9  10 11 12 13 14 15
%     ORDRE DES NOEUDS ICI       : 1 7 2 8 3 9 13 14 15 4  10 5  11 6  12
      SHP(1,1 )=-AUX*DZEM*PAUX2;
      SHP(1,7 )=DEUX*QSI*AUX*DZEM;
      SHP(1,2 )=QSI*DZEM*PAUX;
      SHP(1,8 )=DEUX*QSI*ETA*DZEM;
      SHP(1,3 )=ETA*DZEM*PAUX1;
      SHP(1,9 )=DEUX*ETA*AUX*DZEM;
      SHP(1,13)=AUX*DZEM*DZEP;
      SHP(1,14)=QSI*DZEM*DZEP;
      SHP(1,15)=ETA*DZEM*DZEP;
      SHP(1,4 )=-AUX*DZEP*PAUX5;
      SHP(1,10)=DEUX*QSI*AUX*DZEP;
      SHP(1,5 )=QSI*DZEP*PAUX3;
      SHP(1,11)=DEUX*QSI*ETA*DZEP;
      SHP(1,6 )=ETA*DZEP*PAUX4;
      SHP(1,12)=DEUX*ETA*AUX*DZEP;
    shp1 = SHP(1,:);
    clear SHP;

  otherwise,
    type1
    error('Type of element not implemented yet')
end
