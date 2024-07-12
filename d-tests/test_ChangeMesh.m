clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')

mode1 = 'COPL';
idim = 2;

% On lit les donnees
file1 = 'test_ChangeMesh.inp';
xcrit1 = 1.e-3;
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1, ...
  ListCara1,error1] = ReadMergeAVS2(xcrit1,file1,1);
xcoor1 = xcoor1(:,1:idim);

% Maillage
mail1 = ListMesh1{1};

[mail2,xcoor2] = ChangeMesh2(mail1,'QUAD',xcoor1);
xcoor1 = [xcoor1 ; xcoor2];

file2 = 'resu_ChangeMesh.inp';
ListMesh1{1} = mail2;
error1 = WriteMergeAVS(xcoor1,ListMesh1,ListChpo1,ListnMesh1, ...
                       ListChml1,ListCara1,file2)

disp('VERIFIER CONTENU resu_ChangeMesh.inp')
disp('TEST PASSE AVEC SUCCES')

quit

