% Test de lissage selectif
clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')

mode1 = 'DEPL';
idim = 2;

% On lit les donnees
% """"""""""""""""""
file1 = 'test_SmoothChamno.inp';
nrec1 = 1; xcrit1 = 1.e-6;
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1, ...
  ListCara1,error1] = ReadMergeAVS2(xcrit1,file1,nrec1);
xcoor1 = xcoor1(:,1:idim);

mail1 = ListMesh1{1};
chml1 = ListChml1{1};

clear ListMesh1 ListChpo1 ListnMesh1 ListChml1 ListCara1;

[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);
cham1 = ChmlToCham(chml1,mail1,intg1);
[chamno1,intgno1] = ChamToChamno(cham1,mail1,intg1,xcoor1); 

xcrit1 = 1.e-6; % pas de lissage
chamno2 = SmoothChamno(chamno1,modl1,mail1,intgno1,xcrit1);

xcrit1 = 1.e6; % lissage complet
chamno3 = SmoothChamno(chamno1,modl1,mail1,intgno1,xcrit1);

xcrit1 = 0.09; % lissage partiel
chamno4 = SmoothChamno(chamno1,modl1,mail1,intgno1,xcrit1);

% Verifications
[listComp1,listUnit1] = ListCompCham2(chamno1);
Ncomp1 = length(listComp1);
Nptg1 = Nptg(mail1,intgno1);
numerptg1 = [1:Nptg1*Ncomp1];
mapComp1 = MapCompCham(chamno1,mail1,intgno1,numerptg1,listComp1);
U1 = ChamToVect3(chamno1,mail1,intgno1,numerptg1,mapComp1,listComp1);
U2 = ChamToVect3(chamno2,mail1,intgno1,numerptg1,mapComp1,listComp1);
U3 = ChamToVect3(chamno3,mail1,intgno1,numerptg1,mapComp1,listComp1);
U4 = ChamToVect3(chamno4,mail1,intgno1,numerptg1,mapComp1,listComp1);

nmail1 = ChangeMesh2(mail1,'POI1');
[chpo1,nmail1] = ChamnoToChpo2(chamno1,mail1,intgno1,nmail1,xcoor1,mode1);
[chpo2,nmail2] = ChamnoToChpo2(chamno1,mail1,intgno1,nmail1,xcoor1,mode1, ...
                               'geometric');
[chamno5,intgno5] = ChpoToChamno3(chpo1,nmail1,mail1);
U5 = ChamToVect3(chamno5,mail1,intgno5,numerptg1,mapComp1,listComp1);

err1 = norm(U2 - U1)/norm(U1);
err2 = norm(U3 - U5)/norm(U5);
if (err1 > 1.e-8) | (err2 > 1.e-8)
  err1
  err2
  error('probleme')
end

numer1 = nmail1{1}.MAIL';
[listDdl1,listUnit1] = ListCompChpo2(chpo1);
mapddl1 = MapddlChpo(chpo1,nmail1,numer1,listDdl1);
V1 = ChpoToVect3(chpo1,nmail1,numer1,mapddl1,listDdl1);
V2 = ChpoToVect3(chpo2,nmail2,numer1,mapddl1,listDdl1);
err3 = norm(V2 - V1)/norm(V1);
if (err3 > 1.e-8)
  err3
  error('probleme 2')
end


disp('TEST PASSE AVEC SUCCES')
quit

fid = fopen('test_SmoothChamno.pos','w');
  xcoor2 = [xcoor1 zeros(size(xcoor1,1),1)];
  error1 = Write1GMSH(xcoor2,mail1,chamno2,fid);
fclose(fid);
% gmsh -p test_SmoothChamno.pos
