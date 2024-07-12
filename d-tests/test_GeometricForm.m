% Test the geometric detection of a form from a cloud of nodes
clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')
opt2 = 'CUB8';

function [ino1] = num_no(i,j,k,ni,nj,nk)
  ino1 = i + (j-1)*nj + (k-1)*nj*nk;
end

xcoor = 25.*[0 0 0
             1 0 0
             0 1 0
             1 1 1
             2 1 1
             1 2 1
             2 2 1.5
             0.2 1.5 1.2
             1.5 1.5 1.1];
%xcoor = 25.*[1 1 1
%         2 1 2
%         1 2 1
%         2 2 2
%         0.2 1.5 0.2
%         1.5 1.5 1.5];
%xcoor(:,3) = 0.;
%xcoor(:,3) = 0.02*(xcoor(:,1) .^ 2);
if (1 == 0)
[N,I,G] = GeometricForm(xcoor);
NN = inv(N');
xcoor = xcoor * NN;
xcoor(:,1) = xcoor(:,1) + 50;

xcoor =  25.*[0 0 0
              1 1 0.
              1 2 0.
              2 2 0.
              2 1 0.
              1.5 1.5 0.
              0 2 0
              2 0 0];
xcoor(:,3) = 0.02*(xcoor(:,1) .^ 2);
xcoor(:,3) = 0.02*(xcoor(:,2) .^ 2);
xcoor(:,3) = 0.02*(xcoor(:,2) .* xcoor(:,1));
end

if (1 == 0)
a = 0.05; b = 0.; c = 0.0;
a = 0.0; b = 0.; c = 0.0;
[N,I,G] = GeometricForm(xcoor);
xcoorll = 25*[0. 0. 0.
           0.5 0. 0.
           1. 0. 0.
           1.5 0. 0.
           2.  0. 0.
           0. 0.5 0.
           0.5 0.5 0.
           1. 0.5 0.
           1.5 0.5 0.
           2.  0.5 0.
           0. 1.5 0.
           0.5 1.5 0.
           1. 1.5 0.
           1.5 1.5 0.
           2.  1.5 0.];
z = a*(xcoorll(:,1) .^ 2) + b*(xcoorll(:,2) .^2) + 2.*c*(xcoorll(:,1).*xcoorll(:,2));
xcoorll = [xcoorll(:,1:2) z];
nbno0 = size(xcoorll,1);
N(:,1) = -N(:,1);
det(N)
xcoor = xcoorll * N';
end


[LL,MM,CC,err,N] = GeometricQuadraticFit(xcoor);
disp(['erreur '  num2str(err)])
NN = inv(N');
xcoorl = xcoor * NN;

% Generate the cube for which the zero level is this surface
LnS = [3 6 10 14 25 30 50];
LnF = [3 6 10 14 25 30 50];
LnI = [3 6 10 14 25 30 50];
ni = length(LnS); nj = length(LnF); nk = length(LnI);
lx = LnS; ly = LnF; lz = LnI;
idim = 3;
nbno = ni*nj*nk;
%num_no = inline('i + (j-1)*nj + (k-1)*nj*nk','i','j','k','ni','nj','nk');
xcoor1 = zeros(0,idim);
field2 = zeros(0,1);
for k = 1:nk
  for j = 1:nj
    for i = 1:ni
      ino1 = num_no(i,j,k,ni,nj,nk);
      x = lx(i); y = ly(j); z = lz(k);
      xcoor1(ino1,:) = [x y z];
      field2(ino1,:) = (LL * [x y z]') - ([x y z] * MM * [x y z]') - CC;
    end
  end
end

xcoort1 = [xcoor ; xcoor1];
nbno0 = size(xcoor,1);

switch opt2
  case 'CUB8',
nbel = (ni-1)*(nj-1)*(nk-1);
clear mail1 maile1;
topo1 = zeros(nbel,8);
iel1 = 0;
for k = 1:nk-1
  for j = 1:nj-1
    for i = 1:ni-1
      iel1 = iel1 + 1;
      ino1 = num_no(i,j,k,ni,nj,nk);
      ino2 = num_no(i+1,j,k,ni,nj,nk);
      ino3 = num_no(i+1,j+1,k,ni,nj,nk);
      ino4 = num_no(i,j+1,k,ni,nj,nk);
      ino5 = num_no(i,j,k+1,ni,nj,nk);
      ino6 = num_no(i+1,j,k+1,ni,nj,nk);
      ino7 = num_no(i+1,j+1,k+1,ni,nj,nk);
      ino8 = num_no(i,j+1,k+1,ni,nj,nk);
      topo1(iel1,:) = [ino1 ino2 ino3 ino4 ino5 ino6 ino7 ino8];
    end
  end
end
maile1 = struct('TYPE','CUB8','MAIL',nbno0+topo1);
mail1{1} = maile1;

  case 'TET4',
    topo1 = delaunay3(xcoor1(:,1)',xcoor1(:,2)',xcoor1(:,3)', ...
                      {'Qt', 'QbB', 'Qc', 'Qz'} );
    clear mail1 maile1;
    maile1 = struct('TYPE','TET4','MAIL',nbno0+topo1);
    mail1{1} = maile1;
  otherwise
    opt2
    error('Bad option')
end

[chpo1,nmail1] = VectToChpo2(field2,nbno0+[1:nbno],[1:nbno]',[{'SCAL'}]);

% Local mesh
tri = delaunay(xcoorl(:,1),xcoorl(:,2));
clear mail2;
mail2{1} = struct('TYPE','TRI3','MAIL',tri);

fid = fopen('test_GeometricForm.msh','w');
  ListMesh1{1} = mail2;
  error1 = WriteMeshGMSH(xcoort1,ListMesh1,fid);
fclose(fid);

[chamno1,intg1] = ChpoToChamno3(chpo1,nmail1,mail1);
file1 = 'test_GeometricForm.pos';
fid=fopen(file1,'w');
%  error1 = WriteChamnoGMSH(xcoort1,mail1,chamno1,fid)
  LlistComp1{1} = [{'SCAL'}];
  listnames{1} = 'SCAL';
  error1 = WriteChamnoGMSH2(xcoort1,mail1,chamno1, ...
                            LlistComp1,listnames,fid)
fclose(fid);

disp('TES PASSE AVEC SUCCES')
quit

