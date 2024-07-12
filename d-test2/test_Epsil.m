% Test calcul deformations
clear all
close all

path(path,'../matlabEF') 
path(path,'../matlabEF2') 
path(path,'../matlabUtils') 

mode1 = 'TRID';
idim = 3;

file1 = 'test_Epsil.inp';
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

% Modeles
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);

% % Materiau
% matr2 = ChmlToCham(matr1,mail1,intg1);
% 
% % Rigidite
% [rigi1,bsig1] = Rigi9(modl1,matr2,mail1,intg1,xcoor1,mode1,'GenDualOp');

% Assemblage
numer1 = nmail1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
% K1 = RigiToMatrix(rigi1,modl1,mail1, ...
%                   numer1,mapddlDual1,listDdlDual1, ...
%                   numer1,mapddlPrim1,listDdlPrim1);

% Champ de deplacement
U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrim1,listDdlPrim1);

% % Champ de contraintes
% nbptg1 = Nptg(mail1,intg1);
% numer2 = [1:nbptg1];
% [listComDual1,listComPrim1] = ListCompModl3(modl1);                           
% [mapcomDual1,mapddlDual0] = MapcomB2(modl1,mail1,intg1, ...
%                                      numer2,listComDual1, ...
%                                      numer1,listDdlDual1, ...
%                                      'DUAL');
% BS1 = BToMatrix3(bsig1,modl1,mail1,intg1, ...
%                 numer1,mapddlDual1,listDdlDual1, ...
%                 numer2,mapcomDual1,listComDual1, ...
%                 'DUAL');
% nbddl = size(BS1,1);
% 
% cham1 = ChmlToCham(matr1,mail1,intg1);
% Sig1 = ChamToVect3(cham1,mail1,intg1,numer2,mapcomDual1,listComDual1);
% 
% % Cas classique HPP matriciel
% Fint1 = BS1' * Sig1;

% % Sans matrice, HPP
% numerinv1 = InverseList(numer1,max(numer1));
% Option1 = 'HPP';
% U1 = zeros(nbddl,1);
% Fi1 = Fint(Sig1,mail1,numerinv1,mapddlPrim1, ...
%            intg1,modl1, ...
%            Option1,U1);

numerinv1 = InverseList(numer1,max(numer1));
Option1 = 'HPP';
Epsi1 = Epsil(mail1,xcoor1,numerinv1,mapddlPrim1, ...
              intg1,modl1, ...
              Option1,U1);
Option1 = 'Green-Lagrange';
Epsi2 = Epsil(mail1,xcoor1,numerinv1,mapddlPrim1, ...
              intg1,modl1, ...
              Option1,U1);
Option1 = 'Material-Hencky';
Epsi3 = Epsil(mail1,xcoor1,numerinv1,mapddlPrim1, ...
              intg1,modl1, ...
              Option1,U1);
disp('Verifier les valeurs')

% Verif concordance en HPP
U2 = 1.e-5 * U1;
Option1 = 'HPP';
Epsi1a = Epsil(mail1,xcoor1,numerinv1,mapddlPrim1, ...
              intg1,modl1, ...
              Option1,U2);
Option1 = 'Green-Lagrange';
Epsi2a = Epsil(mail1,xcoor1,numerinv1,mapddlPrim1, ...
              intg1,modl1, ...
              Option1,U2);
Option1 = 'Material-Hencky';
Epsi3a = Epsil(mail1,xcoor1,numerinv1,mapddlPrim1, ...
              intg1,modl1, ...
              Option1,U2);
err2a = max(abs(Epsi2a - Epsi1a))/max(abs(Epsi1a))
err3a = max(abs(Epsi3a - Epsi1a))/max(abs(Epsi1a))


disp('TEST PASSE AVEC SUCCES')
quit
