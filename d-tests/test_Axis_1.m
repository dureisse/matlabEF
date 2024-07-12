% On teste l'axisymetrique en ELASTIQUE, RIGIDITE et MASSE
% et en POREUX, RIGIDITE
% et en FLUIDE, PERMEABILITE et COMPRESSIBILITE
clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')
GlobalVar;
idim = 2;
mode = 'AXIS';

file1 = 'test_Axis_1.inp';
nrec1 = 2; xcrit1 = 1.e-5;
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
  ReadMergeAVS2(xcrit1,file1,nrec1);
xcoor1 = xcoor1(:,1:idim);

mail1 = ListMesh1{1};
chpo1 = ListChpo1{1};
nmail1 = ListnMesh1{1};
matr1 = ListChml1{1};
tmail1 = ListMesh1{2};
tchpo1 = ListChpo1{2};
tnmail1 = ListnMesh1{2};
tmatr1 = ListChml1{2};

% Modeles
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','MASSE',mode,idim);
[modl2,intg2] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode,idim);
[modl3,intg3] = ModlIntg13(mail1,'POREUX','RIGIDITE',mode,idim);
[modl4,intg4] = ModlIntg13(tmail1,'FLUIDE','COMPRESSIBILITE',mode,idim);
[modl5,intg5] = ModlIntg13(tmail1,'FLUIDE','PERMEABILITE',mode,idim);
% Materiaux
mater1 = ChmlToCham(matr1,mail1,intg1);
mater2 = ChmlToCham(matr1,mail1,intg2);
mater3 = ChmlToCham(matr1,mail1,intg3);
mater4 = ChmlToCham(tmatr1,tmail1,intg4);
mater5 = ChmlToCham(tmatr1,tmail1,intg5);
% Matrice de masse
[mass1,bgamma1,bu1,rho1] = Mass7(modl1,mater1,mail1,intg1,xcoor1,mode, ...
                                 'GenDualOp','PrimOp','ConstiOp');
% Matrice de raideur
[rigi1,bsigma1,beps1,hooke1] = Rigi9(modl2,mater2,mail1,intg2,xcoor1,mode, ...
                                 'GenDualOp','PrimOp','ConstiOp');
% Matrice de compressibilite
[comp1,bq1,bp1,cmpr1] = Compress7(modl4,mater4,tmail1,intg4,xcoor1,mode, ...
                                 'GenDualOp','PrimOp','ConstiOp');
% Matrice de raideur couplee rigidite-compressibilite
[rigit1,bsigt1,bepst1,hooket1] = RigiCompress8( ...
                                 modl3,mater3,mail1,intg3,xcoor1,mode, ...
                                 'GenDualOp','PrimOp','ConstiOp');
% Matrice de permeabilite
[perm1,bw1,bz1,hp1] = Perm6(modl5,mater5,tmail1,intg5,xcoor1,mode, ...
                                 'GenDualOp','PrimOp','ConstiOp');

% Assemblages
numer1 = nmail1{1}.MAIL';
tnumer1 = tnmail1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[listDdlDual2,listDdlPrim2] = ListDdlModl2(modl2);
[listDdlDual3,listDdlPrim3] = ListDdlModl2(modl3);
[listDdlDual4,listDdlPrim4] = ListDdlModl2(modl4);
[listDdlDual5,listDdlPrim5] = ListDdlModl2(modl5);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
[mapddlDual2,mapddlPrim2] = MapddlRigi2(modl2,mail1, ...
                                        numer1,listDdlDual2, ...
                                        numer1,listDdlPrim2);
[mapddlDual3,mapddlPrim3] = MapddlRigi2(modl3,mail1, ...
                                        numer1,listDdlDual3, ...
                                        numer1,listDdlPrim3);
[mapddlDual4,mapddlPrim4] = MapddlRigi2(modl4,tmail1, ...
                                        tnumer1,listDdlDual4, ...
                                        tnumer1,listDdlPrim4);
[mapddlDual5,mapddlPrim5] = MapddlRigi2(modl5,tmail1, ...
                                        tnumer1,listDdlDual5, ...
                                        tnumer1,listDdlPrim5);
M1 = RigiToMatrix(mass1,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);
K1 = RigiToMatrix(rigi1,modl2,mail1, ...
                  numer1,mapddlDual2,listDdlDual2, ...
                  numer1,mapddlPrim2,listDdlPrim2);
KT1 = RigiToMatrix(rigit1,modl3,mail1, ...
                  numer1,mapddlDual3,listDdlDual3, ...
                  numer1,mapddlPrim3,listDdlPrim3);

C1 = RigiToMatrix(comp1,modl4,tmail1, ...
                  tnumer1,mapddlDual4,listDdlDual4, ...
                  tnumer1,mapddlPrim4,listDdlPrim4);
H1 = RigiToMatrix(perm1,modl5,tmail1, ...
                  tnumer1,mapddlDual5,listDdlDual5, ...
                  tnumer1,mapddlPrim5,listDdlPrim5);
[listComDual1,listComPrim1] = ListCompModl3(modl1);
[listComDual2,listComPrim2] = ListCompModl3(modl2);
[listComDual3,listComPrim3] = ListCompModl3(modl3);
[listComDual4,listComPrim4] = ListCompModl3(modl4);
[listComDual5,listComPrim5] = ListCompModl3(modl5);
nbptg1 = Nptg(mail1,intg1); numerptg1 = [1:nbptg1];
nbptg2 = Nptg(mail1,intg2); numerptg2 = [1:nbptg2];
nbptg3 = Nptg(mail1,intg3); numerptg3 = [1:nbptg3];
nbptg4 = Nptg(tmail1,intg4); numerptg4 = [1:nbptg4];
nbptg5 = Nptg(tmail1,intg5); numerptg5 = [1:nbptg5];
[mapcomDual1,mapddlDual0] = MapcomB2(modl1,mail1,intg1, ...
                                     numerptg1,listComDual1, ...
                                     numer1,listDdlDual1, ...
                                     'DUAL');
[mapcomDual2,mapddlDual0] = MapcomB2(modl2,mail1,intg2, ...
                                     numerptg2,listComDual2, ...
                                     numer1,listDdlDual2, ...
                                     'DUAL');
[mapcomDual3,mapddlDual0] = MapcomB2(modl3,mail1,intg3, ...
                                     numerptg3,listComDual3, ...
                                     numer1,listDdlDual3, ...
                                     'DUAL');
[mapcomDual4,mapddlDual0] = MapcomB2(modl4,tmail1,intg4, ...
                                     numerptg4,listComDual4, ...
                                     tnumer1,listDdlDual4, ...
                                     'DUAL');
[mapcomDual5,mapddlDual0] = MapcomB2(modl5,tmail1,intg5, ...
                                     numerptg5,listComDual5, ...
                                     tnumer1,listDdlDual5, ...
                                     'DUAL');
[mapcomPrim1,mapddlPrim0] = MapcomB2(modl1,mail1,intg1, ...
                                     numerptg1,listComPrim1, ...
                                     numer1,listDdlPrim1, ...
                                     'PRIMAL');
[mapcomPrim2,mapddlPrim0] = MapcomB2(modl2,mail1,intg2, ...
                                     numerptg2,listComPrim2, ...
                                     numer1,listDdlPrim2, ...
                                     'PRIMAL');
[mapcomPrim3,mapddlPrim0] = MapcomB2(modl3,mail1,intg3, ...
                                     numerptg3,listComPrim3, ...
                                     numer1,listDdlPrim3, ...
                                     'PRIMAL');
[mapcomPrim4,mapddlPrim0] = MapcomB2(modl4,tmail1,intg4, ...
                                     numerptg4,listComPrim4, ...
                                     tnumer1,listDdlPrim4, ...
                                     'PRIMAL');
[mapcomPrim5,mapddlPrim0] = MapcomB2(modl5,tmail1,intg5, ...
                                     numerptg5,listComPrim5, ...
                                     tnumer1,listDdlPrim5, ...
                                     'PRIMAL');

BG1 = BToMatrix3(bgamma1,modl1,mail1,intg1, ...
                 numer1,mapddlDual1,listDdlDual1, ...
                 numerptg1,mapcomDual1,listComDual1, ...
                 'DUAL');
BU1 = BToMatrix3(bu1,modl1,mail1,intg1, ...
                 numer1,mapddlPrim1,listDdlPrim1, ...
                 numerptg1,mapcomPrim1,listComPrim1, ...
                 'PRIMAL');
RH1 = DToMatrix2(rho1,modl1,mail1,intg1, ...
                 numerptg1,mapcomPrim1,listComPrim1, ...
                 numerptg1,mapcomDual1,listComDual1);

BS1 = BToMatrix3(bsigma1,modl2,mail1,intg2, ...
                 numer1,mapddlDual2,listDdlDual2, ...
                 numerptg2,mapcomDual2,listComDual2, ...
                 'DUAL');
BE1 = BToMatrix3(beps1,modl2,mail1,intg2, ...
                 numer1,mapddlPrim2,listDdlPrim2, ...
                 numerptg2,mapcomPrim2,listComPrim2, ...
                 'PRIMAL');
HO1 = DToMatrix2(hooke1,modl2,mail1,intg2, ...
                 numerptg2,mapcomPrim2,listComPrim2, ...
                 numerptg2,mapcomDual2,listComDual2);

BQ1 = BToMatrix3(bq1,modl4,tmail1,intg4, ...
                 tnumer1,mapddlDual4,listDdlDual4, ...
                 numerptg1,mapcomDual4,listComDual4, ...
                 'DUAL');
BP1 = BToMatrix3(bp1,modl4,tmail1,intg4, ...
                 tnumer1,mapddlPrim4,listDdlPrim4, ...
                 numerptg4,mapcomPrim4,listComPrim4, ...
                 'PRIMAL');
CO1 = DToMatrix2(cmpr1,modl4,tmail1,intg4, ...
                 numerptg4,mapcomPrim4,listComPrim4, ...
                 numerptg4,mapcomDual4,listComDual4);

BST1 = BToMatrix3(bsigt1,modl3,mail1,intg3, ...
                 numer1,mapddlDual3,listDdlDual3, ...
                 numerptg3,mapcomDual3,listComDual3, ...
                 'DUAL');
BET1 = BToMatrix3(bepst1,modl3,mail1,intg3, ...
                 numer1,mapddlPrim3,listDdlPrim3, ...
                 numerptg3,mapcomPrim3,listComPrim3, ...
                 'PRIMAL');
HOT1 = DToMatrix2(hooket1,modl3,mail1,intg3, ...
                 numerptg3,mapcomPrim3,listComPrim3, ...
                 numerptg3,mapcomDual3,listComDual3);

BW1 = BToMatrix3(bw1,modl5,tmail1,intg5, ...
                 tnumer1,mapddlDual5,listDdlDual5, ...
                 numerptg5,mapcomDual5,listComDual5, ...
                 'DUAL');
BZ1 = BToMatrix3(bz1,modl5,tmail1,intg5, ...
                 tnumer1,mapddlPrim5,listDdlPrim5, ...
                 numerptg5,mapcomPrim5,listComPrim5, ...
                 'PRIMAL');
HP1 = DToMatrix2(hp1,modl5,tmail1,intg5, ...
                 numerptg5,mapcomPrim5,listComPrim5, ...
                 numerptg5,mapcomDual5,listComDual5);

% Verifications of consistency
% """"""""""""""""""""""""""""
% Mass from associated operators
M2 = BG1' * RH1 * BU1;
err1 = norm(full(M1 - M2)) / norm(full(M1));
if (err1 > 1e-7)
  err1
  error('erreur 1')
end
% Stiffness from associated operators
K2 = BS1' * HO1 * BE1;
err2 = norm(full(K1 - K2)) / norm(full(K1));
if (err2 > 1e-7)
  err2
  error('erreur 2')
end
% Stiffness from particular assembling of coupled problem
K2 = RigiToMatrix(rigit1,modl3,mail1, ...
                  numer1,mapddlDual2,listDdlDual2, ...
                  numer1,mapddlPrim2,listDdlPrim2);
err2a = norm(full(K1 - K2)) / norm(full(K1));
if (err2a > 1e-7)
  err2a
  error('erreur 2a')
end
% Compressibility from associated operators
C2 = BQ1' * CO1 * BP1;
err3 = norm(full(C1 - C2)) / norm(full(C1));
if (err3 > 1e-7)
  err3
  error('erreur 3')
end
% Compressibility from particular assembling of coupled problem
C2 = RigiToMatrix(rigit1,modl3,mail1, ...
                  tnumer1,mapddlDual4,listDdlDual4, ...
                  tnumer1,mapddlPrim4,listDdlPrim4);
C2 = -1. * C2;
err3a = norm(full(C1 - C2)) / norm(full(C1));
if (err3a > 1e-7)
  err3a
  error('erreur 3a')
end
% Permeability from associated operators
H2 = BW1' * HP1 * BZ1;
err4 = norm(full(H1 - H2)) / norm(full(H1));
if (err1 > 1e-7)
  err1
  error('erreur 4')
end

% partie diagonale seule
ddlp = mapddlPrim3(:,3)';
ddlp = ddlp(find(ddlp));
nddl = size(KT1,2);
ddlu = setdiff([1:nddl],ddlp);
KT2 = BST1' * HOT1 * BET1;
err3 = norm(full(KT1(ddlu,ddlu) - KT2(ddlu,ddlu))) / norm(full(KT1(ddlu,ddlu)));
err3 = err3 + ...
       norm(full(KT1(ddlp,ddlp) - KT2(ddlp,ddlp))) / norm(full(KT1(ddlp,ddlp)));
if (err3 > 1e-7)
  err3
  error('erreur 3')
end

% Verifications of values
% """""""""""""""""""""""
U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrim1,listDdlPrim1);
XVAL1 = U1' * M1 * U1;
if (abs(XVAL1 - 8.4410442861)/8.4410442861) > 1.E-6
  XVAL1
  error('mauvaise valeur 1')
end
XVAL2 = U1' * K1 * U1;
if (abs(XVAL2 - 19.703142489)/19.703142489) > 1.E-6
  XVAL2
  error('mauvaise valeur 2')
end
U3 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrim3,listDdlPrim3);
XVAL3 = U3' * KT1 * U3;
if (abs(XVAL3 + 0.66152630382)/ 0.66152630382) > 1.E-6
  XVAL3
  error('mauvaise valeur 3')
end
U4 = ChpoToVect3(tchpo1,tnmail1,tnumer1,mapddlPrim5,listDdlPrim5);
XVAL4 = U4' * H1 * U4;
if (abs(XVAL4 - 2.5617873746)/2.5617873746) > 1.E-6
  XVAL4
  error('mauvaise valeur 4')
end

disp('TEST PASSE AVEC SUCCES')
quit
