% On teste l'axisymetrique en ELASTIQUE, RIGIDITE et MASSE
% en 1D plonge en 2D
clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')
GlobalVar;
idim = 2;
mode = 'AXIS';

file1 = 'test_Axis_2.inp';
nrec1 = 1; xcrit1 = 1.e-5;
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
  ReadMergeAVS2(xcrit1,file1,nrec1);
xcoor1 = xcoor1(:,1:idim);

mail1 = ListMesh1{1};
chpo1 = ListChpo1{1};
nmail1 = ListnMesh1{1};
matr1 = ListChml1{1};

% Modeles
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','MASSE',mode,idim);
% Materiaux
mater1 = ChmlToCham(matr1,mail1,intg1);
% Matrice de masse
[mass1,bgamma1,bu1,rho1] = Mass7(modl1,mater1,mail1,intg1,xcoor1,mode, ...
                                 'GenDualOp','PrimOp','ConstiOp');

% Assemblages
numer1 = nmail1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
%[listDdlDual2,listDdlPrim2] = ListDdlModl2(modl2);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
%[mapddlDual2,mapddlPrim2] = MapddlRigi2(modl2,mail1, ...
%                                        numer1,listDdlDual2, ...
%                                        numer1,listDdlPrim2);
M1 = RigiToMatrix(mass1,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);
[listComDual1,listComPrim1] = ListCompModl3(modl1);
nbptg1 = Nptg(mail1,intg1); numerptg1 = [1:nbptg1];
[mapcomDual1,mapddlDual0] = MapcomB2(modl1,mail1,intg1, ...
                                     numerptg1,listComDual1, ...
                                     numer1,listDdlDual1, ...
                                     'DUAL');
[mapcomPrim1,mapddlPrim0] = MapcomB2(modl1,mail1,intg1, ...
                                     numerptg1,listComPrim1, ...
                                     numer1,listDdlPrim1, ...
                                     'PRIMAL');

BG1 = BToMatrix3(bgamma1,modl1,mail1,intg1, ...
                 numer1,mapddlDual1,listDdlDual1, ...
                 numerptg1,mapcomDual1,listComDual1, ...
                 'DUAL');
BU1 = BToMatrix3(bu1,modl1,mail1,intg1, ...
                 numer1,mapddlPrim1,listDdlPrim1, ...
                 numerptg1,mapcomPrim1,listComPrim1, ...
                 'PRIMAL');
RH1 = DToMatrix2(rho1,modl1,mail1,intg1, ...
                 numerptg1,mapcomPrim1,listComPrim1, ...
                 numerptg1,mapcomDual1,listComDual1);


% Verifications of consistency
% """"""""""""""""""""""""""""
% Mass from associated operators
M2 = BG1' * RH1 * BU1;
err1 = norm(full(M1 - M2)) / norm(full(M1));
if (err1 > 1e-7)
  err1
  error('erreur 1')
end

% Verifications of values
% """""""""""""""""""""""
ALPHA = 30. * pi / 180.;
F1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlDual1,listDdlDual1);
P1 = zeros(size(M1,2),1);
  for no = 1:length(numer1)
    no1 = numer1(no);
    ddl1 = mapddlDual1(no,1);
    z1 = xcoor1(no1,2);
    P1(ddl1,:) = z1 * cos(ALPHA);
    P1(ddl1+1,:) = z1 * sin(ALPHA);
  end
F2 = M1 * P1;
err4 = norm(F1 - F2) / norm(F1);
if (err4 > 1e-7)
  err4
  error('erreur 4')
end

disp('TEST PASSE AVEC SUCCES')
quit
