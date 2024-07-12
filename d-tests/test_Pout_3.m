% Test poutre d'Euler-Bernoulli - courbure negligee
% raccord de poutres
clear all
close all
path(path,'../matlabEF') 
path(path,'../matlabUtils') 

mode1 = 'POUT';
idim = 2;

file1 = 'test_Pout_3.inp';
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

% Modele
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);

% Materiau
matr2 = ChmlToCham(matr1,mail1,intg1);

% On cree le champ par elements aux noeuds des tangentes
[chamno1,intgno1] = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1);

%% On le lisse selon un critere de brisure de tangentes
%xcrit1 = 1.e-2;
%SmoothChamno(chamno1,modl1,mail1,intgno1,xcrit1)

[rigi2,bsigma2,bepsilon2,d2] = Rigi9(modl1,matr2,mail1,intg1,xcoor1,mode1,...
                                     'GenDualOp','PrimOp','ConstiOp');
[rigi3,bsigma3,bepsilon3,d3] = Rigi9(modl1,matr2,mail1,intg1,xcoor1,mode1,...
                                     chamno1,'GenDualOp','PrimOp','ConstiOp');
[r1,modlr1] = IntegrOperator2(modl1,mail1,intg1,xcoor1,mode1);

% Assemblage
numer1 = RenumberRcm(mail1,xcoor1);
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[listComDual1,listComPrim1] = ListCompModl3(modl1);
[listComDualr1,listComPrimr1] = ListCompModl3(modlr1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
nbptg1 = Nptg(mail1,intg1);
mapcomPrimr1 = [1:nbptg1]';
mapcomDualr1 = mapcomPrimr1;
numer2 = [1:nbptg1];
[mapcomDual1,mapddlDual0] = MapcomB2(modl1,mail1,intg1, ...
                                     numer2,listComDual1, ...
                                     numer1,listDdlDual1, ...
                                     'DUAL');
[mapcomPrim1,mapddlPrim0] = MapcomB2(modl1,mail1,intg1, ...
                                     numer2,listComPrim1, ...
                                     numer1,listDdlPrim1, ...
                                     'PRIMAL');

K1 = RigiToMatrix(rigi2,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);
B1 = BToMatrix3(bsigma2,modl1,mail1,intg1, ...
                numer1,mapddlDual1,listDdlDual1, ...
                numer2,mapcomDual1,listComDual1, ...
                'DUAL');
bb1 = BToMatrix3(bepsilon2,modl1,mail1,intg1, ...
                 numer1,mapddlPrim1,listDdlPrim1, ...
                 numer2,mapcomPrim1,listComPrim1, ...
                 'PRIMAL');
dd1 = DToMatrix2(d2,modl1,mail1,intg1, ...
                 numer2,mapcomPrim1,listComPrim1, ...
                 numer2,mapcomDual1,listComDual1);
K3 = RigiToMatrix(rigi3,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);
B3 = BToMatrix3(bsigma3,modl1,mail1,intg1, ...
                numer1,mapddlDual1,listDdlDual1, ...
                numer2,mapcomDual1,listComDual1, ...
                'DUAL');
bb3 = BToMatrix3(bepsilon3,modl1,mail1,intg1, ...
                 numer1,mapddlPrim1,listDdlPrim1, ...
                 numer2,mapcomPrim1,listComPrim1, ...
                 'PRIMAL');
dd3 = DToMatrix2(d3,modl1,mail1,intg1, ...
                 numer2,mapcomPrim1,listComPrim1, ...
                 numer2,mapcomDual1,listComDual1);
R1 = DToMatrix2(r1,modlr1,mail1,intg1,...
                numer2,mapcomPrimr1,listComPrimr1, ...
                numer2,mapcomDualr1,listComDualr1);

K2 = B1' * dd1 * bb1;

% Champ par point
U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrim1,listDdlPrim1);

% Verification
XVAL1 = U1' * K1 * U1;
if (abs(XVAL1 - 1.3349325825)/1.3349325825 > 1.E-10)
  XVAL1
  error('mauvaise valeur')
end

err1 = norm(full(K1 - K2)) / norm(full(K1));
if (err1 > 1e-7)
  err1
  error('erreur')
end

err2 = norm(full(K1 - K3)) / norm(full(K1))
err3 = norm(full(B1 - B3)) / norm(full(B1))
err4 = norm(full(bb1 - bb3)) / norm(full(bb1))
err5 = norm(full(dd1 - dd3)) / norm(full(dd1))

disp('TEST PASSE AVEC SUCCES')
quit
