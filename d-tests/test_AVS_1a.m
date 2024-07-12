clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')

% On lit les donnees
disp('Read AVS-UCD format')
file1 = 'test_AVS_1b.inp';
xcrit1 = -1.e-6;
disp('Attention, pb quand on merge les noeuds sur le champ par elements')
disp('Pour voir le pb : passer xcrit1 a 1.e-6')
[xcoort1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
  ReadMergeAVS2(xcrit1,file1,1);

% On ecrit les resultats
disp('Write AVS-UCD format')
file2 = 'test_AVS_1ba.inp';
error1 = WriteMergeAVS(xcoort1,ListMesh1,ListChpo1,ListnMesh1, ...
                       ListChml1,ListCara1,file2);

disp('FAIRE UN diff test_AVS_1b.inp test_AVS_1ba.inp')
system('diff -b test_AVS_1b.inp test_AVS_1ba.inp')

disp('Write GMSH-Mesh 1.0 format')
fid = fopen('test_GMSH_1b.msh','w');
  error1 = WriteMeshGMSH(xcoort1,ListMesh1,fid);
fclose(fid);

disp('Write GMSH-Postprocessing scalar and vector point file')
mail1 = ListMesh1{1};
chpo0 = ListChpo1{1};
nmail0 = ListnMesh1{1};
[chamno1,junk] = ChpoToChamno3(chpo0,nmail0,mail1);     
fid = fopen('test_GMSH_1b.pos','w');
%  error1 = WriteChamnoGMSH(xcoort1,mail1,chamno1,fid)
  LlistComp1{1} = [{'SCAL'}];
  listnames{1} = 'SCAL';
  error1 = WriteChamnoGMSH2(xcoort1,mail1,chamno1, ...
                            LlistComp1,listnames,fid)
fclose(fid);
% ! gmsh -p test_GMSH_1b.pos

quit


[chamno2,junk] = ChpoToChamno3(chpo0,nmail0,mail1);
clear Lchamno1; Lchamno1{1} = chamno2;
nmail2 = nmail0;
clear Lchpo2; Lchpo2{1} = chpo0;
listComp2 = [{'SCAL'}];
fid = fopen('test_GMSH_1c.pos','w');
  [error1] = WriteLChGMSH3(fid,'test',xcoort1,[0.], ...
                            mail1,Lchamno1,listComp2, ...
                            nmail2,Lchpo2,listComp2)
fclose(fid);
% ! gmsh -p test_GMSH_1c.pos

