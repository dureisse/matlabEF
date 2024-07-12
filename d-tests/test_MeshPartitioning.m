clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')
%path(path,'meshpart')
 
idim = 2;
   
% On lit les donnees
file1 = 'test_AVS_1.inp';
xcrit1 = 1.e-6;
[xcoort1,ListMesh1,ListChpo1,ListnMesh1,ListChml1, ...
  ListCara1,error1] = ReadMergeAVS2(xcrit1,file1,1);
  xcoort1 = xcoort1(:,1:idim);

mail1 = ListMesh1{1};

%n = 2;
%cham1 = MeshPartitioning(mail1,n,xcoort1)

lnsdm1 = [3 2];
theta = 0.5;
Q = [cos(theta) -sin(theta); sin(theta) cos(theta)];
cham1 = MeshPartitioning0(mail1,lnsdm1,Q,xcoort1);

clear ListMesh2 ListChpo2 ListnMesh2 ListChml2 ListCara2;
file3 = 'test_MeshPartitioning.inp';
ListMesh2{1} = mail1;
ListChpo2{1} = [];
ListnMesh2{1} = [];
ListChml2{1} = cham1;
ListCara2{1} = [];
error1 = WriteMergeAVS(xcoort1,ListMesh2,ListChpo2,ListnMesh2, ...
                       ListChml2,ListCara2,file3);

quit
