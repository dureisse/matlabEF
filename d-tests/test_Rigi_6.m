% Test de l'element CUB8 en 3D
clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')
GlobalVar;
mode1 = 'TRID';
idim = 3;

% On lit les donnees
file1 = 'test_Rigi_6.inp';
xcrit1 = -1.;
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1, ...
  ListCara1,error1] = ReadMergeAVS2(xcrit1,file1,1);
xcoor1 = xcoor1(:,1:idim);

% Maillage
mail1 = ListMesh1{1};

% Champ de deplacement
chpo1 = ListChpo1{1};
nmail1 = ListnMesh1{1};

% Modele
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);

nbelt1 = Nelements2(mail1);
nbptg1 = Nptg(mail1,intg1);
numer2 = [1:nbptg1];

% Materiau
mater2 = ListChml1{1};
mater2 = ChmlToCham(mater2,mail1,intg1);

% Rigidite et Bsigma et B
[rigi2,d2,bsig2,b2] = Rigi9(modl1,mater2,mail1,intg1,xcoor1,mode1, ...
                            'ConstiOp','GenDualOp','PrimOp');
[r1,modlr1] = IntegrOperator2(modl1,mail1,intg1,xcoor1,mode1);

numer1 = nmail1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[listComDual1,listComPrim1] = ListCompModl3(modl1);                           
[listComDualr1,listComPrimr1] = ListCompModl3(modlr1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
[mapcomDual1,mapddlDual0] = MapcomB2(modl1,mail1,intg1, ...
                                     numer2,listComDual1, ...
                                     numer1,listDdlDual1, ...
                                     'DUAL');
[mapcomPrim1,mapddlPrim0] = MapcomB2(modl1,mail1,intg1, ...
                                     numer2,listComPrim1, ...
                                     numer1,listDdlPrim1, ...
                                     'PRIMAL');
mapcomPrimr1 = [1:nbptg1]';
mapcomDualr1 = mapcomPrimr1;

U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrim1,listDdlPrim1);
Fsig1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlDual1,listDdlDual1);
K1 = RigiToMatrix(rigi2,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);

B1 = BToMatrix3(bsig2,modl1,mail1,intg1, ...
                numer1,mapddlDual1,listDdlDual1, ...
                numer2,mapcomDual1,listComDual1, ...
                'DUAL');
bb1 = BToMatrix3(b2,modl1,mail1,intg1, ...
                 numer1,mapddlPrim1,listDdlPrim1, ...
                 numer2,mapcomPrim1,listComPrim1, ...
                 'PRIMAL');
dd1 = DToMatrix2(d2,modl1,mail1,intg1, ...
                 numer2,mapcomPrim1,listComPrim1, ...
                 numer2,mapcomDual1,listComDual1);
R1 = DToMatrix2(r1,modlr1,mail1,intg1,...
                numer2,mapcomPrimr1,listComPrimr1, ...
                numer2,mapcomDualr1,listComDualr1);

K2 = B1' * dd1 * bb1;
err1 = norm(full(K1 - K2)) / norm(full(K1));

if (err1 > 1e-7)
  err1
  error('erreur')
end

eneref1 =  477.02529455;

ener1 = U1' * K1 * U1;
err2 = abs(ener1 - eneref1) / eneref1;
if (err2 > 1e-7)
  err2
  error('error 2')
end

epsi1 = bb1 * U1;
sigm1 = dd1 * epsi1;
F1 = B1'*sigm1;
F2 = K1 * U1;
err3 = norm(F1 - F2) / norm(F1);
if (err3 > 1e-7)
  err3
  error('error 3')
end

Epsi1 = epsi1(mapcomPrim1);
Sigm1 = sigm1(mapcomDual1);
Ener1 = sum(Sigm1 .* Epsi1,2);
ener2 = ones(size(Ener1))' * R1 * Ener1;
err4 = abs(ener2 - eneref1) / eneref1;
if (err4 > 1e-7)
  err4
  error('error 4')
end

disp('TEST PASSE AVEC SUCCES')
quit