clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')
%mode1 = 'COPL'; % Plane Stress
%mode1 = 'DEPL'; % Plane Strain
mode1 = 'AXIS'; % Axisymetrique
idim = 2;

% Recuperation des donnees
% """"""""""""""""""""""""
disp('Recuperation des donnees')
fid = fopen('test_Poreux_4.inp','rt');
  [xcoor1,mail1,chpo1,chml1,cara1,error1] = Read1AVS6(fid);
fclose(fid);
xcoor1 = xcoor1(:,1:idim);

nbno1 = size(xcoor1,1);
clear nmail1;
  nmail1{1} = struct('MAIL',[1:nbno1]','TYPE','POI1');
numer1 = nmail1{1}.MAIL';


% Modele couple
% """""""""""""
[modl1,intg1] = ModlIntg13(mail1,'POREUX','RIGIDITE',mode1,idim);

nbelt1 = Nelements2(mail1);
nbptg1 = Nptg(mail1,intg1);
numerptg1 = [1:nbptg1];

% Coefficients materiau
% """""""""""""""""""""
mater1 = ChmlToCham(chml1,mail1,intg1);

% Matrice de rigidite/compressibilite couplee
% """""""""""""""""""""""""""""""""""""""""""
disp('Matrice de rigidite/compressibilite couplee')
[rigit1,bsigt1,bepst1] = RigiCompress8(modl1,mater1,mail1,intg1, ...
                                       xcoor1,mode1,'GenDualOp','PrimOp');

% Assemblage (sur tous les ddl : modl1)
% Vont donc ensembles :
% (numer1,mapddlDualt1,listDdlDualt1)
% (numer1,mapddlPrimt1,listDdlPrimt1)
% (numerptg1,mapcomDualt1,listComDualt1)
% (numerptg1,mapcomPrimt1,listComPrimt1)
disp('Assemblage')
[listDdlDualt1,listDdlPrimt1] = ListDdlModl2(modl1);
[listComDualt1,listComPrimt1] = ListCompModl3(modl1);
[mapddlDualt1,mapddlPrimt1] = MapddlRigi2(modl1,mail1, ...
                                          numer1,listDdlDualt1, ...
                                          numer1,listDdlPrimt1);
[mapcomDualt1,mapddlDualt0] = MapcomB2(modl1,mail1,intg1, ...
                                       numerptg1,listComDualt1, ...
                                       numer1,listDdlDualt1, ...
                                       'DUAL');
[mapcomPrimt1,mapddlPrimt0] = MapcomB2(modl1,mail1,intg1, ...
                                       numerptg1,listComPrimt1, ...
                                       numer1,listDdlPrimt1, ...
                                       'PRIMAL');
KT1 = RigiToMatrix(rigit1,modl1,mail1, ...
                   numer1,mapddlDualt1,listDdlDualt1, ...
                   numer1,mapddlPrimt1,listDdlPrimt1);
BSIGT1 = BToMatrix3(bsigt1,modl1,mail1,intg1, ...
                    numer1,mapddlDualt1,listDdlDualt1, ...
                    numerptg1,mapcomDualt1,listComDualt1, ...
                    'DUAL');
BEPST1 = BToMatrix3(bepst1,modl1,mail1,intg1, ...
                    numer1,mapddlPrimt1,listDdlPrimt1, ...
                    numerptg1,mapcomPrimt1,listComPrimt1, ...
                    'PRIMAL');


% Modele de permeabilite
% """"""""""""""""""""""
mail2 = ChangeMesh2(mail1,'TRI3');
[modl2,intg2] = ModlIntg13(mail2,'FLUIDE','PERMEABILITE',mode1,idim);

% Matrice de permeabilite
% """""""""""""""""""""""
disp('Matrice de permeabilite')
[perm2,bsig2,b2] = Perm6(modl2,mater1,mail2,intg2,xcoor1,mode1, ...
                         'GenDualOp','PrimOp');

% Assemblage (sur tous les ddl : modl1)
disp('Assemblage')
HH1 = RigiToMatrix(perm2,modl2,mail2, ...
                   numer1,mapddlDualt1,listDdlDualt1, ...
                   numer1,mapddlPrimt1,listDdlPrimt1);


% Champ par point
U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrimt1,listDdlPrimt1); 
XVAL1 = U1' * KT1 * U1;
%if (abs(XVAL1 - 5.5040416008)/5.5040416008 > 1.E-9)
if (abs(XVAL1 +  0.20382324028)/ 0.20382324028 > 1.E-9)
  XVAL1
  error('mauvaise valeur')
end
XVAL2 = U1' * HH1 * U1;
%if (abs(XVAL2 - 251.8518222222)/251.8518222222 > 1.E-8)
if (abs(XVAL2 -  3552.3976242249229)/ 3552.3976242249229 > 1.E-8)
  XVAL2
  error('mauvaise valeur')
end



disp('TEST PASSE AVEC SUCCES')
quit
