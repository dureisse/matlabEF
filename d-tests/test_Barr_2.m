clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')

mode1 = 'BARR';
idim  = 2;

% On lit les donnees
% """"""""""""""""""
file1 = 'test_Barr_2.inp';
nrec1 = 5; xcrit1 = 1.e-6;
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1, ...
  ListCara1,error1] = ReadMergeAVS2(xcrit1,file1,nrec1);
xcoor1 = xcoor1(:,1:idim);

% Maillage SEG2 des poutres
mail1 = ListMesh1{1};
% Materiau
mater1 = ListChml1{1};

% Efforts exterieurs
cfext1 = ListChpo1{2};
nfext1 = ListnMesh1{2};

% Deplacements imposes UY
cdepi1 = ListChpo1{3};
ndepi1 = ListnMesh1{3};

% Deplacements imposes UX
cdepi2 = ListChpo1{4};
ndepi2 = ListnMesh1{4};

% Solution de reference
cref1 = ListChpo1{5};
nref1 = ListnMesh1{5};

% Modele
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);

% Rigidite
mater2 = ChmlToCham(mater1,mail1,intg1);
[cara2,junk]  = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1);
[rigi1,bsigma1,bepsilon1] = Rigi9(modl1,mater2,mail1,intg1,xcoor1,mode1, ...
                                  cara2,'GenDualOp','PrimOp');

% Assemblage
% """"""""""
numer1 = RenumberRcm(mail1,xcoor1);
% nmail1 = ChangeMesh2(mail1,'POI1');
% numer1 = nmail1{1}.MAIL';
% clear nmail1;

[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
K1 = RigiToMatrix(rigi1,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);

U1 = ChpoToVect3(cref1,nref1,numer1,mapddlPrim1,listDdlPrim1);
XVAL1 = U1' * K1 * U1;
if (abs(XVAL1 - 7.1038220528)/7.1038220528 > 1.E-8)
  XVAL1
  error('mauvaise valeur')
end

disp('TEST PASSE AVEC SUCCES')
quit

