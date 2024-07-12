clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')

disp('Read AVS-UCD format')
fid = fopen('test_AVS_1.inp','rt');
  [xcoor1,mail1,chpo1,chml1,cara1,error1] = Read1AVS6(fid);
fclose(fid);

disp('Write AVS-UCD format')
fid = fopen('test_AVS_1a.inp','w');
  error1 = Write1AVS5(xcoor1,mail1,chpo1,chml1,cara1,fid);
fclose(fid);

disp('FAIRE UN diff test_AVS_1.inp test_AVS_1a.inp')
system('diff -b test_AVS_1.inp test_AVS_1a.inp')

disp('Write GMSH-Mesh 1.0 format')
fid = fopen('test_GMSH_1.msh','w');
  ListMesh1{1} = mail1;
  error1 = WriteMeshGMSH(xcoor1,ListMesh1,fid);
fclose(fid);

quit
