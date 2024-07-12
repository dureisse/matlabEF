% On teste la matrice de masse sur un SEG2 en 3D
clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')
GlobalVar;
idim = 3;
mode = 'BARR';

file1 = 'test_Mass_4.inp';
nrec1 = 1; xcrit1 = 1.e-5;
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
  ReadMergeAVS2(xcrit1,file1,nrec1);
xcoor1 = xcoor1(:,1:idim);

mail1 = ListMesh1{1};
chpo1 = ListChpo1{1};
nmail1 = ListnMesh1{1};
matr1 = ListChml1{1};

% Modele
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','MASSE',mode,idim);

% Materiau
mater1 = ChmlToCham(matr1,mail1,intg1);

% Matrice de masse
[mass1,bgamma1,bu1,rho1] = Mass7(modl1,mater1,mail1,intg1,xcoor1,mode, ...
                                 'GenDualOp','PrimOp','ConstiOp');

% Assemblages
numer1 = nmail1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
M1 = RigiToMatrix(mass1,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);

% Verification
U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrim1,listDdlPrim1);                 
XVAL1 = U1' * M1 * U1;
if (abs(XVAL1 - 1774.0514531105003)/1774.0514531105003 > 1.E-6)
  XVAL1
  error('mauvaise valeur')
end                

disp('TEST PASSE AVEC SUCCES')
quit
