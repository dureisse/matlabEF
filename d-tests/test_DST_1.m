clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')

mode1 = 'DSHE';
idim = 3;

% On lit les donnees
file1 = 'test_DST_1.inp';
xcrit1 = -1.;
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1, ...
  ListCara1,error1] = ReadMergeAVS2(xcrit1,file1,1);
xcoor1 = xcoor1(:,1:idim);

% Maillage
mail1 = ListMesh1{1};
 
% Champ de deplacement
chpo1 = ListChpo1{1};
nmail1 = ListnMesh1{1};

% Materiau
matr1 = ListChml1{1};

% Modele
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);

% Materiau
matr2 = ChmlToCham(matr1,mail1,intg1);

% Base locale
[chamno1,intgno1] = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1);

% Rigidites elementaires
rigi1 = Rigi9(modl1,matr2,mail1,intg1,xcoor1,mode1,chamno1,'Spin',1.e-7);

% Assemblages
numer1 = nmail1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);

U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrim1,listDdlPrim1);
K1 = RigiToMatrix(rigi1,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);

% Verification
XVAL1 = U1' * K1 * U1;
% Elem de ref, E = 1, nu = 0.3, epai = 1.
err1 = abs(XVAL1 - 1.7658262514632226);
% Elem de ref, E = 1, nu = 0.3, epai = 0.1
err1 = abs(XVAL1 - 0.1436462260812054);
% Elem de ref, E = 0.5, nu = 0., epai = 0.2
err1 = abs(XVAL1 - 0.1027667275385838);
% Elem de ref x 4, E = 1, nu = 0.3, epai = 0.1
err1 = abs(XVAL1 - 17.0625394735513);
% Elem distordu, E = 1, nu = 0.3, epai = 0.1
err1 = abs(XVAL1 - 1.007360763118792);
if (err1 > 1.e-8)
  error('probleme')
end

disp('TEST PASSE AVEC SUCCES')
quit
