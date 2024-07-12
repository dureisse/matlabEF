clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')

% On lit les donnees
disp('Read AVS-UCD format')
file1 = 'test_AVS_2.inp';
xcrit1 = 1.e-6;
%[xcoort1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
%  ReadMergeAVS2(xcrit1,file1,1);
[xcoort1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
  ReadMergeAVS2(xcrit1,file1,2);

% On ecrit les resultats
disp('Write AVS-UCD format')
file2 = 'test_AVS_2a.inp';
error1 = WriteMergeAVS(xcoort1,ListMesh1,ListChpo1,ListnMesh1, ...
                       ListChml1,ListCara1,file2);

disp('FAIRE UN diff test_AVS_2.inp test_AVS_2a.inp')
system('diff -b test_AVS_2.inp test_AVS_2a.inp')

disp('Write GMSH-Mesh 1.0 format') 
fid = fopen('test_GMSH_2a.msh','w');
  error1 = WriteMeshGMSH(xcoort1,ListMesh1,fid);
fclose(fid);

disp('Write GMSH-Postprocessing scalar and vector point file')
mail1 = ListMesh1{1};
chpo0 = ListChpo1{1};
nmail0 = ListnMesh1{1};
[chamno1,junk] = ChpoToChamno3(chpo0,nmail0,mail1);     
%nmail1 = ListnMesh1{2};
%chpo1  = ListChpo1{2};
name1  = 'test_AVS_2a';
ListScalarComp1 = [{'UX'}];
ListVectorComp1 = [{'UX'} {'UY'} {'UZ'}];
fid = fopen('test_GMSH_2a.pos','w');
%  error1 = WriteChpoGMSH(xcoort1,nmail1,chpo1,fid,name1, ...
%                         ListScalarComp1,ListVectorComp1)
% Champ de scalaire aux points seul
%  error1 = WriteChpoGMSH(xcoort1,nmail1,chpo1,fid,name1, ...
%                         ListScalarComp1)
% Champ de scalaire par elements
%  error1 = WriteChamnoGMSH(xcoort1,mail1,chamno1,fid)  
  LlistComp1{1} = [{'SCAL'}];
  listnames{1} = 'SCAL';
  error1 = WriteChamnoGMSH2(xcoort1,mail1,chamno1, ...
                            LlistComp1,listnames,fid)
% Champ de vecteur aux points seul
%  error1 = WriteChpoGMSH(xcoort1,nmail1,chpo1,fid,name1, ...
%                         ListVectorComp1)
%  error1 = WriteChGMSH(xcoort1,mail1,chamno1, ...
%                       nmail1,chpo1,fid,name1, ...
%                       ListVectorComp1)
fclose(fid);
% ! gmsh -p test_GMSH_2a.pos

quit
