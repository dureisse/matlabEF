% teste le passage chpo vers chamno en discret et en direct
clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')

idim = 2;

fid = fopen('test_Vect_1.inp','rt');
  [xcoor1,mail1,chpo1,chml1,cara1,error1] = Read1AVS6(fid);
fclose(fid);
xcoor1 = xcoor1(:,1:idim);

% En direct
nbno1 = size(xcoor1,1);
clear nmail1;
  nmail1{1} = struct('MAIL',[1:nbno1]','TYPE','POI1');
[chamno1,intg1] = ChpoToChamno3(chpo1,nmail1,mail1);
nptg1 = Nptg(mail1,intg1);

% En matriciel
numer1 = [1:nbno1];
numerptg = [1:nptg1];
T1 = ChpoToChamnoOperator(numer1,mail1,numerptg);

% On commence par passer en matriciel
[listDdl1,listUnit1] = ListCompChpo2(chpo1);
mapddl1 = MapddlChpo(chpo1,nmail1,numer1,listDdl1);
U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddl1,listDdl1);

[listComp1,listUnit0] = ListCompCham2(chamno1);
mapComp1 = MapCompCham(chamno1,mail1,intg1,numerptg,listComp1);    
C1 = ChamToVect3(chamno1,mail1,intg1,numerptg,mapComp1,listComp1);

% On verifie
U2 = U1(mapddl1);
C2 = T1 * U2;
C3 = zeros(size(C1));
C3(mapComp1) = C2;

err1 = norm(C1 - C3) / norm(C1);
if (err1 > 1.e-7)
  err1
  error('erreur')
end

disp('TEST PASSE AVEC SUCCES')
quit
