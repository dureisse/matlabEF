clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')

fid = fopen('test_Vect_1.inp','rt');
  [xcoor1,mail1,chpo1,chml1,cara1,error1] = Read1AVS6(fid);
fclose(fid);

xcoor1 = xcoor1(:,1:2);
nbno1 = size(xcoor1,1);
clear nmail1;
  nmail1{1} = struct('MAIL',[1:nbno1]','TYPE','POI1');

numer1 = [1:nbno1];
[listComp1,listUnit1] = ListCompChpo2(chpo1);
mapddl1 = MapddlChpo(chpo1,nmail1,numer1,listComp1);
U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddl1,listComp1);

[chpo2,nmail2] = VectToChpo2(U1,numer1,mapddl1,listComp1, ...
                             listUnit1); 

err1 = norm(chpo1{1}{1}.XVAL - chpo2{1}{1}.XVAL) / ...
	norm(chpo1{1}{1}.XVAL);

if (err1 > 1.e-7)
  err1
  error('erreur')
end

disp('TEST PASSE AVEC SUCCES')
quit
