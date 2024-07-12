% On teste l'intersection de maillages

clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')

idim = 2;

% On lit les donnees
file1 = 'test_AVS_2.inp';
xcrit1 = 1.e-6;
[xcoort1,ListMesh1,ListChpo1,ListnMesh1,ListChml1, ...
  ListCara1,error1] = ReadMergeAVS2(xcrit1,file1,2);
xcoort1 = xcoort1(:,1:idim);

mail1 = ListMesh1{1};
mail2 = ListMesh1{2};

% mail3 devrait etre mail1
mail3 = IntersectMesh(mail1,mail1);
[lnod1,lind1] = MeshToGraph(mail1);
[lnod3,lind3] = MeshToGraph(mail3);
err1 = norm(lnod1-lnod3)/norm(lnod1);
err2 = norm(lind1-lind3)/norm(lind1);
if (err1 > 1.e-10) | (err2 > 1.e-10)
  err1
  err2
  error('PB')
end

% mail3 devrait etre mail1
mail3 = IntersectMesh(mail1,mail2);
[lnod1,lind1] = MeshToGraph(mail1);
[lnod3,lind3] = MeshToGraph(mail3);
err1 = norm(lnod1-lnod3)/norm(lnod1);
err2 = norm(lind1-lind3)/norm(lind1);
if (err1 > 1.e-10) | (err2 > 1.e-10)
  err1
  err2
  error('PB2')
end

disp('TEST PASSE AVEC SUCCES')
quit
