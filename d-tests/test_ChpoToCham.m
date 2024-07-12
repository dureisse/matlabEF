% Test du changement de champ par point en champ par element
% (interpolation)
clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')

mode1 = 'COPL';
idim = 2;

file1 = 'test_Cham.inp';
xcrit1 = 1.e-6;
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
  ReadMergeAVS2(xcrit1,file1,1);
xcoor1 = xcoor1(:,[1:idim]);

% Champ par point defini sur nmail1,mail1
chpo1 = ListChpo1{1};
mail1 = ListMesh1{1};
nmail1 = ListnMesh1{1};

% Modele
[modl1,intg1] = ModlIntg13(mail1,'FLUIDE','COMPRESSIBILITE',mode1,idim);

[cham1] = ChpoToCham(chpo1,nmail1,mail1,intg1);

% Version matricielle
numer1 = nmail1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
					                    numer1,listDdlPrim1);

[Chpo1] = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrim1,listDdlPrim1);

Nptg1 = Nptg(mail1,intg1);
numerptg1 = [1:Nptg1];
[T1] = ChpoToChamOperator(numer1,mail1,numerptg1,intg1);

Cham1 = T1 * Chpo1;

[listComDual1,listComPrim1] = ListCompModl3(modl1);
[mapcomPrim1,mapddlPrim0] = MapcomB2(modl1,mail1,intg1, ...
                                     numerptg1,listComPrim1, ...
				                     numer1,listDdlPrim1, ...
				                     'PRIMAL');

[Cham2] = ChamToVect3(cham1,mail1,intg1,numerptg1,mapcomPrim1,listComPrim1);

err1 = sqrt(norm(Cham1-Cham2)/norm(Cham2));
if (err1 > 1.e-8)
  err1
  error('PB')
end

disp('TEST PASSE AVEC SUCCES')
quit
