% On teste les QUA8
clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')
idim = 2;

fid = fopen('test_Mass_1b.inp','rt');
  [xcoor1,mail1,chpo1,chml1,cara1,error1] = Read1AVS6(fid);
fclose(fid);
xcoor1 = xcoor1(:,1:idim);

nmail1 = ChangeMesh2(mail1,'POI1');

% Modele
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','MASSE','COPL',idim);

% Materiau
mater1 = ManuChml(mail1,'RHO','',100.);
mater1 = ChmlToCham(mater1,mail1,intg1);

% Matrice de masse
%[mass1,bgamma1,bu1] = Mass6(modl1,mater1,mail1,intg1,xcoor1,'COPL');
[mass1,bgamma1,bu1] = Mass7(modl1,mater1,mail1,intg1,xcoor1,'COPL', ...
                            'GenDualOp','PrimOp');

mater2 = ManuChml(mail1,'YOUN','',1.,'RHO','',100.,'NU','',0.3);
mater2 = ChmlToCham(mater2,mail1,intg1);

%[mass2,bgamma2,bu2] = Mass6(modl1,mater2,mail1,intg1,xcoor1,'COPL');
[mass2,bgamma2,bu2] = Mass7(modl1,mater2,mail1,intg1,xcoor1,'COPL', ...
                            'GenDualOp','PrimOp');

err1 = max(max(max(abs(mass1{1}.XVAL-mass2{1}.XVAL)))) / ...
       max(max(max(abs(mass1{1}.XVAL))))

numer1 = nmail1{1}.MAIL';
[listCompDual1,listCompPrimal1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrimal1] = MapddlRigi2(modl1,mail1, ...
                                          numer1,listCompDual1, ...
                                          numer1,listCompPrimal1);
M1 = RigiToMatrix(mass2,modl1,mail1, ...
		    numer1,mapddlDual1,listCompDual1, ...
		    numer1,mapddlPrimal1,listCompPrimal1);

% Masse
V = ones(size(M1,1),1);
mm1 = V' * M1 * V

% Pseudo-inertie 2D
chpo2 = CoorMail(nmail1,xcoor1);
[listComp1,listUnit1] = ListCompChpo2(chpo2);
mapddl1 = MapddlChpo(chpo2,nmail1,numer1,listComp1);
U1 = ChpoToVect3(chpo2,nmail1,numer1,mapddl1,listComp1);
nbddl1 = size(U1,1);
ix = findoccur({'UX'},listComp1);
iy = findoccur({'UY'},listComp1);
UX1 = U1(mapddl1(:,ix),1);
UY1 = U1(mapddl1(:,iy),1);
V1 = zeros(nbddl1,2);
V1(mapddl1(:,ix),1) = UY1;
V1(mapddl1(:,ix),2) = -1. * UX1;
inert1 = V1' * M1 * V1

% Champ par point
U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrimal1,listCompPrimal1); 
XVAL1 = U1' * M1 * U1;
%if (abs(XVAL1 - 213.33333)/213.33333 > 1.E-8)
if (abs(XVAL1 - 218.666666)/218.666666 > 1.E-8)
  XVAL1
  error('mauvaise valeur')
end

disp('TEST PASSE AVEC SUCCES')
quit
