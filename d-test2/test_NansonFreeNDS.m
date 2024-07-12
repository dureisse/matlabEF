% Test 
clear all
close all

path(path,'../matlabEF') 
path(path,'../matlabEF2') 
path(path,'../matlabUtils') 

mode1 = 'TRID';
idim = 3;

file1 = 'test_NansonFreeNDS.inp';
nrec1 = 1; % nombre d'enregistrements
xcrit1 = -1.e-6; % critere de proximite de noeuds
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
  ReadMergeAVS2(xcrit1,file1,nrec1);
xcoort1 = xcoor1(:,1:idim);
mail1s = ListMesh1{1};
chpo1s = ListChpo1{1};
nmail1s = ListnMesh1{1};
clear ListMesh1 ListChpo1 ListnMesh1 ListChml1 ListCara1;

% Integration nodale ! sur l'element omega_0
intgno1s = SegmentIntgNo(mail1s,'nodes');
% Modele
[modl1s,intg1s] = ModlIntg13(mail1s,'ELASTIQUE','MASSE',mode1,idim,intgno1s);

% Pour assemblage
numer1s = nmail1s{1}.MAIL';
nmax1 = max(numer1s);
numerinv1s = InverseList(numer1s,nmax1);
% Vecteur position des noeuds
X1s = xcoor1(numer1s,:);

% Vecteur deplacement
[listDdlDual1s,listDdlPrim1s] = ListDdlModl2(modl1s);
[mapddlDual1s,mapddlPrim1s] = MapddlRigi2(modl1s,mail1s, ...
                                          numer1s,listDdlDual1s, ...
                                          numer1s,listDdlPrim1s);
U1 = ChpoToVect3(chpo1s,nmail1s,numer1s,mapddlPrim1s,listDdlPrim1s);
U1s = U1(mapddlPrim1s);
% Vecteur forces generalisees
F1 = ChpoToVect3(chpo1s,nmail1s,numer1s,mapddlDual1s,listDdlDual1s);

% Configuration courante
% X1s = X1s + U1s;


tic
[nds1] = NansonFreeNDS(mail1s,modl1s,intg1s, ...
                       X1s,numerinv1s);
toc

NDS1 = 0*U1;
NDS1(mapddlPrim1s) = nds1;
norm(F1-NDS1)/norm(F1)

[chpo2s,nmail2s] = VectToChpo2(NDS1,numer1s,mapddlPrim1s,listDdlPrim1s);
fid = fopen('test_NansonFreeNDS.pos','w');
  chamno1s = [];
  LlistComp1s = [];
  clear LlistComp2; LlistComp2{1} = [{'UX'} {'UY'} {'UZ'}];
  error1 = WriteChGMSH2(fid,'Normales',xcoort1, ...
                                mail1s,chamno1s,LlistComp1s, ...
                                nmail1s,chpo2s,LlistComp2)
fclose(fid);
fid = fopen('test_NansonFreeNDSm.pos','w');
clear ListMesh1; ListMesh1{1} = mail1s;
  error1 = WriteMeshGMSH(xcoort1,ListMesh1,fid)
fclose(fid);

disp('TEST PASSE AVEC SUCCES')
quit
