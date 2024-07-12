% On teste la matrice de compressibilite sur le bord
clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')
GlobalVar;
mode1 = 'COPL';
idim = 2;

fid = fopen('test_Compress_3.inp','rt');
  [xcoor1,mail1,chpo1,chml1,cara1,error1] = Read1AVS6(fid);
fclose(fid);
xcoor1 = xcoor1(:,1:idim);

nbno1 = size(xcoor1,1);
clear nmail1;
  nmail1{1} = struct('MAIL',[1:nbno1]','TYPE','POI1');

% Modele
[modl1,intg1] = ModlIntg13(mail1,'FLUIDE','COMPRESSIBILITE',mode1,idim);

% Materiau
mater1 = ManuChml(mail1,'MOB','',1.);
mater1 = ChmlToCham(mater1,mail1,intg1);

% Matrice de compressibilite bord
compres2 = Compress7(modl1,mater1,mail1,intg1,xcoor1,mode1);

% Assemblage
numer1 = nmail1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
S1 = RigiToMatrix(compres2,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);

% Verification
lthermique = [{'Q'}];
F1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlDual1,lthermique);

flux1 = ManuChpo(nmail1,'P','',1.);
P1 = ChpoToVect3(flux1,nmail1,numer1,mapddlPrim1,listDdlPrim1);
F2 = S1 * P1;
err1 = norm(F1 - (F2 / 1.))/ norm(F1);
if (err1 > 1.e-7)
  error('erreur')
end

disp('TEST PASSE AVEC SUCCES')
quit
