% Flambement de poutre
clear all
close all

path(path,'../matlabEF') 
path(path,'../matlabEF2') 
path(path,'../matlabUtils') 

mode1 = 'POUT';
idim = 3;

file1 = 'test_Pout3D_6.inp';
nrec1 = 1; % nombre d'enregistrements
xcrit1 = -1.e-6; % critere de proximite de noeuds
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
  ReadMergeAVS2(xcrit1,file1,nrec1);
xcoor1 = xcoor1(:,1:idim);

mail1 = ListMesh1{1};
chpo1 = ListChpo1{1};
nmail1 = ListnMesh1{1};
matr1 = ListChml1{1};

clear ListMesh1 ListChpo1 ListnMesh1 ListChml1 ListCara1;

% Modeles
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);

% % On cree le champ par elements aux noeuds des tangentes
% % mail1 est ecrase par son reoriente eventuel
% xcoor2 = [0. -1.];
% [chamno1,intgno1,mail1] = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1,xcoor2);

% Materiau
matr2 = ChmlToCham(matr1,mail1,intg1);

rigi1 = Rigi9(modl1,matr2,mail1,intg1,xcoor1,mode1);

% Assemblage
% numer1 = RenumberRcm(mail1,xcoor1);
numer1 = nmail1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
K1 = RigiToMatrix(rigi1,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);

% Deplacements
U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrim1,listDdlPrim1);

% Verificat1Gion
XVAL1 = U1' * K1 * U1;
if (abs(XVAL1 - 5287.3298077)/5287.329807 > 1.E-6)
  XVAL1
  error('mauvaise valeur')
end


disp('TEST PASSE AVEC SUCCES')
quit
