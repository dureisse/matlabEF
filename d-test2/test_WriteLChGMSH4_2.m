clear all;
close all;
path(path,'../matlabEF')
path(path,'../matlabEF2')
path(path,'../matlabUtils')
GlobalVar;

file1 = 'test_WriteLChGMSH4_2.inp';
nrec1 = 1; xcrit1 = 1.e-5;
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
  ReadMergeAVS2(xcrit1,file1,nrec1);
mail1 = ListMesh1{1};
chpo1 = ListChpo1{1};
nmail1 = ListnMesh1{1};
chml1 = ListChml1{1};

% On cree des champs par element définis aux noeuds
[chamno1,intg1] = ChpoToChamno3(chpo1,nmail1,mail1);
chamno2 = ChmlToChamno(chml1,mail1);
chamno1{1}{1}.XVAL = chamno1{1}{1}.XVAL + chamno2{1}{1}.XVAL;
chamno1{1}{2}.XVAL = chamno1{1}{2}.XVAL + chamno2{1}{2}.XVAL;
chamno1{1}{3}.XVAL = chamno1{1}{3}.XVAL + chamno2{1}{3}.XVAL;
chamno1{1}{4}.XVAL = chamno1{1}{4}.XVAL + chamno2{1}{4}.XVAL;

% chamno1{1}{4}.XVAL = 0 * chamno1{1}{3}.XVAL; % une 4e composante
% chamno1{1}{4}.COMP = 'ZERO';
% chamno1{1}{4}.UNIT = '';
% 
% chpo1{1}{4}.XVAL = 0 * chpo1{1}{3}.XVAL; % une 4e composante
% chpo1{1}{4}.COMP = 'ZERO';
% chpo1{1}{4}.UNIT = '';

% On cree les listes de champs(t)
lt1 = [0:1:3];
nt1 = length(lt1);
clear Lchpo1 Lchamno1;
for it1 = 1:nt1
    Lchpo1{it1} = chpo1;
    Lchamno1{it1} = chamno1;
    for i1 = 1:4
        xval1 = chpo1{1}{i1}.XVAL;
        yval1 = chamno1{1}{i1}.XVAL;
        xval2 = (1 + lt1(it1)*0.1)*xval1;
        yval2 = (1 + lt1(it1)*0.1)*yval1;
        Lchpo1{it1}{1}{i1}.XVAL = xval2;
        Lchamno1{it1}{1}{i1}.XVAL = yval2;
    end
end

file2 = 'test_WriteLChGMSH4_2a.msh';
fid = fopen(file2,'w');
  name1 = 'test';
  name2 = 'testch';
  listComp2 = [{'UX'}]; % scalaire
  error1 = WriteLChGMSH4(fid,xcoor1,lt1, ...
                         mail1,Lchamno1,listComp2,name2, ...
                         nmail1,Lchpo1,listComp2,name1)
fclose(fid)


file2 = 'test_WriteLChGMSH4_2.msh';
fid = fopen(file2,'w');
  name1 = 'test';
  name2 = 'testch';
  listComp2 = [{'UX'} {'UY'} {'UZ'}]; % vecteur
  error1 = WriteLChGMSH4(fid,xcoor1,lt1, ...
                         mail1,Lchamno1,listComp2,name2, ...
                         nmail1,Lchpo1,listComp2,name1)
fclose(fid)

file2 = 'test_WriteLChGMSH4_2b.msh';
fid = fopen(file2,'w');
  name1 = 'test';
  name2 = 'testch';
  listComp2 = [{'UY'} {'ZERO'} {'ZERO'} ...
               {'ZERO'} {'UY'} {'ZERO'} ...
               {'ZERO'} {'ZERO'} {'UY'}]; % tenseur
  error1 = WriteLChGMSH4(fid,xcoor1,lt1, ...
                         mail1,Lchamno1,listComp2,name2, ...
                         nmail1,Lchpo1,listComp2,name1)
  listComp1 = listComp2;
  nmail2 = nmail1;
  Lchpo2 = Lchpo1;
fclose(fid)

disp('TEST PASSE AVEC SUCCES')
quit