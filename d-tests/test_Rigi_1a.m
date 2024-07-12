% Idem test_Rigi_1 mais en deformations planes
clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')
mode1 = 'DEPL';
idim = 2;

fid = fopen('test_Rigi_1.inp','rt');
  [xcoor1,mail1,chpo1,chml1,cara1,error1] = Read1AVS6(fid);
fclose(fid);
xcoor1 = xcoor1(:,1:idim);

nbno1 = size(xcoor1,1);
clear nmail1;
  nmail1{1} = struct('MAIL',[1:nbno1]','TYPE','POI1');

% Modele
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);
nbelt1 = size(mail1{1}.MAIL,1);

% Materiau
clear mater2;
mater2{1} = struct('COMP','YOUN','UNIT','','XVAL',1.*ones(nbelt1,1));
mater2{2} = struct('COMP','RHO','UNIT','','XVAL',100.*ones(nbelt1,1));
mater2{3} = struct('COMP','NU','UNIT','','XVAL',0.3*ones(nbelt1,1));
% mater2 = ManuChml(mail1,'YOUN','',1.,'RHO','',100.,'NU','',0.3);
mater2 = ChmlToCham(mater2,mail1,intg1);

% Rigidite
[rigi2,bsigma2,bepsilon2] = Rigi9(modl1,mater2,mail1,intg1,xcoor1,mode1,...
                                  'GenDualOp','PrimOp');

numer1 = nmail1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                          numer1,listDdlDual1, ...
                                          numer1,listDdlPrim1);
K1 = RigiToMatrix(rigi2,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);

% Champ par point
U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrim1,listDdlPrim1); 
XVAL1 = U1' * K1 * U1;
if (abs(XVAL1 - 3.8461538462)/3.8461538462 > 1.E-10)
  XVAL1
  error('mauvaise valeur')
end

disp('TEST PASSE AVEC SUCCES')
quit
