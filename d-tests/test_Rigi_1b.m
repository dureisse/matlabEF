% On teste les QUA8
clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')
idim = 2;

fid = fopen('test_Rigi_1b.inp','rt');
  [xcoor1,mail1,chpo1,chml1,cara1,error1] = Read1AVS6(fid);
fclose(fid);
xcoor1 = xcoor1(:,1:idim);

nmail1 = ChangeMesh2(mail1,'POI1');

% Modele
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE','COPL',idim);

% Materiau
mater2 = ManuChml(mail1,'YOUN','',1.,'RHO','',100.,'NU','',0.3);
mater2 = ChmlToCham(mater2,mail1,intg1);

% Matrice de rigidite
[rigi2,bsigma2,beps2] = Rigi9(modl1,mater2,mail1,intg1,xcoor1,'COPL', ...
                            'GenDualOp','PrimOp');

numer1 = nmail1{1}.MAIL';
[listCompDual1,listCompPrimal1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrimal1] = MapddlRigi2(modl1,mail1, ...
                                          numer1,listCompDual1, ...
                                          numer1,listCompPrimal1);
K1 = RigiToMatrix(rigi2,modl1,mail1, ...
		    numer1,mapddlDual1,listCompDual1, ...
		    numer1,mapddlPrimal1,listCompPrimal1);

% Champ par point
U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrimal1,listCompPrimal1); 
XVAL1 = U1' * K1 * U1;
if (abs(XVAL1 - 10.256410256)/10.256410256 > 1.E-8)
  XVAL1
  error('mauvaise valeur')
end

disp('TEST PASSE AVEC SUCCES')
quit
