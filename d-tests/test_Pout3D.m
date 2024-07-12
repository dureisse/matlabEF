% Test poutre d'Euler-Bernoulli - courbure negligee : rigidite et masse
clear all
close all
path(path,'../matlabEF') 
path(path,'../matlabUtils') 

mode1 = 'POUT';
idim = 3;

file1 = 'test_Pout3D.inp';
nrec1 = 1; % nombre d'enregistrements
xcrit1 = -1.e-6; % critere de proximite de noeuds
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
  ReadMergeAVS2(xcrit1,file1,nrec1);
xcoor1 = xcoor1(:,1:idim);

mail1 = ListMesh1{1};
chpo1 = ListChpo1{1};
matr1 = ListChml1{1};
nmail1 = ChangeMesh2(mail1,'POI1');

% Modele
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);
if (1 == 0)
[modl2,intg2] = ModlIntg13(mail1,'ELASTIQUE','MASSE',mode1,idim);
end

% Materiau
matr2 = ChmlToCham(matr1,mail1,intg1);
if (1 == 0)
matr3 = ChmlToCham(matr1,mail1,intg2);
end

[rigi2,bsigma2,bepsilon2,d2] = Rigi9(modl1,matr2,mail1,intg1,xcoor1,mode1,...
                                     'GenDualOp','PrimOp','ConstiOp');
if (1 == 0)
[mass2,bgamma2,bu2,rho2] = Mass7(modl2,matr3,mail1,intg2,xcoor1,mode1,...
                                 'GenDualOp','PrimOp','ConstiOp');
end
[r1,modlr1] = IntegrOperator2(modl1,mail1,intg1,xcoor1,mode1);
if (1 == 0)
[r2,modlr2] = IntegrOperator2(modl2,mail1,intg2,xcoor1,mode1);
end

% Assemblage
nmail1 = ChangeMesh2(mail1,'POI1');
numer1 = nmail1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[listComDual1,listComPrim1] = ListCompModl3(modl1);
if (1 == 0)
[listComDual2,listComPrim2] = ListCompModl3(modl2);
end
[listComDualr1,listComPrimr1] = ListCompModl3(modlr1);
if (1 == 0)
[listComDualr2,listComPrimr2] = ListCompModl3(modlr2);
end
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
nbptg1 = Nptg(mail1,intg1);
if (1 == 0)
nbptg2 = Nptg(mail1,intg2);
end
mapcomPrimr1 = [1:nbptg1]';
mapcomDualr1 = mapcomPrimr1;
if (1 == 0)
mapcomPrimr2 = [1:nbptg2]';
mapcomDualr2 = mapcomPrimr2;
end
numerptg1 = [1:nbptg1];
if (1 == 0)
numerptg2 = [1:nbptg2];
end
[mapcomDual1,mapddlDual0] = MapcomB2(modl1,mail1,intg1, ...
                                     numerptg1,listComDual1, ...
                                     numer1,listDdlDual1, ...
                                     'DUAL');
[mapcomPrim1,mapddlPrim0] = MapcomB2(modl1,mail1,intg1, ...
                                     numerptg1,listComPrim1, ...
                                     numer1,listDdlPrim1, ...
                                     'PRIMAL');
if (1 == 0)
[mapcomDual2,mapddlDual0] = MapcomB2(modl2,mail1,intg2, ...
                                     numerptg2,listComDual2, ...
                                     numer1,listDdlDual1, ...
                                     'DUAL');
[mapcomPrim2,mapddlPrim0] = MapcomB2(modl2,mail1,intg2, ...
                                     numerptg2,listComPrim2, ...
                                     numer1,listDdlPrim1, ...
                                     'PRIMAL');
end

K1 = RigiToMatrix(rigi2,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);
B1 = BToMatrix3(bsigma2,modl1,mail1,intg1, ...
                numer1,mapddlDual1,listDdlDual1, ...
                numerptg1,mapcomDual1,listComDual1, ...
                'DUAL');
bb1 = BToMatrix3(bepsilon2,modl1,mail1,intg1, ...
                 numer1,mapddlPrim1,listDdlPrim1, ...
                 numerptg1,mapcomPrim1,listComPrim1, ...
                 'PRIMAL');
dd1 = DToMatrix2(d2,modl1,mail1,intg1, ...
                 numerptg1,mapcomPrim1,listComPrim1, ...
                 numerptg1,mapcomDual1,listComDual1);
if (1 == 0)
M2 = RigiToMatrix(mass2,modl2,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);
Bg2 = BToMatrix3(bgamma2,modl2,mail1,intg2, ...
                numer1,mapddlDual1,listDdlDual1, ...
                numerptg2,mapcomDual2,listComDual2, ...
                'DUAL');
Bu2 = BToMatrix3(bu2,modl2,mail1,intg2, ...
                 numer1,mapddlPrim1,listDdlPrim1, ...
                 numerptg2,mapcomPrim2,listComPrim2, ...
                 'PRIMAL');
Rho2 = DToMatrix2(rho2,modl2,mail1,intg2, ...
                 numerptg2,mapcomPrim2,listComPrim2, ...
                 numerptg2,mapcomDual2,listComDual2);
end
R1 = DToMatrix2(r1,modlr1,mail1,intg1,...
                numerptg1,mapcomPrimr1,listComPrimr1, ...
                numerptg1,mapcomDualr1,listComDualr1);
if (1 == 0)
R2 = DToMatrix2(r2,modlr2,mail1,intg2,...
                numerptg2,mapcomPrimr2,listComPrimr2, ...
                numerptg2,mapcomDualr2,listComDualr2);
end

K2 = B1' * dd1 * bb1;
if (1 == 0)
M1 = Bg2' * Rho2 * Bu2;
end

% Champ par point
U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrim1,listDdlPrim1);

% Verification
XVAL1 = U1' * K1 * U1;
if (abs(XVAL1 - 147.77996026)/147.77996026 > 1.E-10)
  XVAL1
  error('mauvaise valeur')
end

err1 = norm(full(K1 - K2)) / norm(full(K1));
if (err1 > 1e-7)
  err1
  error('erreur')
end

if (1 == 0)
XVAL2 = U1' * M2 * U1
% Impossible de comparer a Cast3M qui calcule la matrice de masse
% de facon analytique et pas numerique !
%if (abs(XVAL2 - 93.8094292804)/93.8094292804 > 1.E-10)
%  XVAL2
%  error('mauvaise valeur 2')
%end

err2 =  norm(full(M1 - M2)) / norm(full(M2));
if (err1 > 1e-7)
  err1
  error('erreur 2')
end
end


disp('TEST PASSE AVEC SUCCES')
quit
