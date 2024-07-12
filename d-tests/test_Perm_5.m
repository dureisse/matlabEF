% Test de l'element TET4 en 3D
clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')
GlobalVar;
mode1 = 'TRID';
idim = 3;

% On lit les donnees
file1 = 'test_Perm_5.inp';
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
[modl1,intg1] = ModlIntg13(mail1,'FLUIDE','PERMEABILITE',mode1,idim);

nbelt1 = Nelements2(mail1);
nbptg1 = Nptg(mail1,intg1);
numer2 = [1:nbptg1];

% Materiau isotrope
mater2 = ListChml1{1};
mater2 = ChmlToCham(mater2,mail1,intg1);

% Materiau faussement anisotrope
PERM3 = ExtrCham2(mater2,{'PERM'}); PERM3 = PERM3{1}{1}.XVAL;
VISC3 = ExtrCham2(mater2,{'VISC'}); VISC3 = VISC3{1}{1}.XVAL;
clear mater3 matere3;
matere3{1} = struct('COMP','PEXX','UNIT','','XVAL',PERM3);
matere3{2} = struct('COMP','PEYY','UNIT','','XVAL',PERM3);
matere3{3} = struct('COMP','PEZZ','UNIT','','XVAL',PERM3);
matere3{4} = struct('COMP','PEXY','UNIT','','XVAL',0.*PERM3);
matere3{5} = struct('COMP','PEYX','UNIT','','XVAL',0.*PERM3);
matere3{6} = struct('COMP','PEYZ','UNIT','','XVAL',0.*PERM3);
matere3{7} = struct('COMP','PEZY','UNIT','','XVAL',0.*PERM3);
matere3{8} = struct('COMP','PEZX','UNIT','','XVAL',0.*PERM3);
matere3{9} = struct('COMP','PEXZ','UNIT','','XVAL',0.*PERM3);
matere3{10} = struct('COMP','VISC','UNIT','','XVAL',VISC3);
mater3{1} = matere3;
clear matere3 PERM3 VISC3;

% Permeabilite et associes
[rigi2,d2,bsig2,b2] = Perm6(modl1,mater2,mail1,intg1,xcoor1,mode1, ...
                            'ConstiOp','GenDualOp','PrimOp');
[rigi3,d3,bsig3,b3] = Perm6(modl1,mater3,mail1,intg1,xcoor1,mode1, ...
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
K2 = RigiToMatrix(rigi3,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);

B1 = BToMatrix3(bsig2,modl1,mail1,intg1, ...
                numer1,mapddlDual1,listDdlDual1, ...
                numer2,mapcomDual1,listComDual1, ...
                'DUAL');
B2 = BToMatrix3(bsig3,modl1,mail1,intg1, ...
                numer1,mapddlDual1,listDdlDual1, ...
                numer2,mapcomDual1,listComDual1, ...
                'DUAL');
bb1 = BToMatrix3(b2,modl1,mail1,intg1, ...
                 numer1,mapddlPrim1,listDdlPrim1, ...
                 numer2,mapcomPrim1,listComPrim1, ...
                 'PRIMAL');
bb2 = BToMatrix3(b3,modl1,mail1,intg1, ...
                 numer1,mapddlPrim1,listDdlPrim1, ...
                 numer2,mapcomPrim1,listComPrim1, ...
                 'PRIMAL');
dd1 = DToMatrix2(d2,modl1,mail1,intg1, ...
                 numer2,mapcomPrim1,listComPrim1, ...
                 numer2,mapcomDual1,listComDual1);
dd2 = DToMatrix2(d3,modl1,mail1,intg1, ...
                 numer2,mapcomPrim1,listComPrim1, ...
                 numer2,mapcomDual1,listComDual1);
R1 = DToMatrix2(r1,modlr1,mail1,intg1,...
                numer2,mapcomPrimr1,listComPrimr1, ...
                numer2,mapcomDualr1,listComDualr1);

err1a = norm(full(K1 - K2)) / norm(full(K1));
if (err1a > 1e-7)
  err1a
  error('erreur')
end
err2a = norm(full(B1 - B2)) / norm(full(B1));
if (err2a > 1e-7)
  err2a
  error('erreur')
end
err3a = norm(full(bb1 - bb2)) / norm(full(bb1));
if (err3a > 1e-7)
  err3a
  error('erreur')
end
err4a = norm(full(dd1 - dd2)) / norm(full(dd1));
if (err4a > 1e-7)
  err4a
  error('erreur')
end

K2 = B1' * dd1 * bb1;
err1 = norm(full(K1 - K2)) / norm(full(K1));
if (err1 > 1e-7)
  err1
  error('erreur')
end

eneref1 = 835.94790923249645;

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

ener3 = Sigm1' * R1 * Epsi1;
ener3 = ener3(1,1) + ener3(2,2) + ener3(3,3);
err5 = abs(ener3 - eneref1) / eneref1;
if (err5 > 1e-7)
  err5
  error('error 5')
end


disp('TEST PASSE AVEC SUCCES')
quit
