clear all
close all

path(path,'../matlabEF')
path(path,'../matlabUtils')

mode1 = 'DEPL'; % Deformations planes
idim = 2;
 
% Recuperation des donnees
% """"""""""""""""""""""""
disp('====> donnees')
file1 = 'test_Assembl.inp';
nrec1 = 1; xcrit1 = 1.e-5;
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
  ReadMergeAVS2(xcrit1,file1,nrec1);
xcoor1 = xcoor1(:,1:idim);

%   Maillage
mail1 = ListMesh1{1};
chpo1 = ListChpo1{1};
nmail1 = ListnMesh1{1};
chml1 = ListChml1{1};

% Modele
% """"""
disp('====> modele couple')
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);

% Materiau
% """"""""
mater1 = ChmlToCham(chml1,mail1,intg1);

% Rigidite
% """"""""
rigit1 = Rigi9(modl1,mater1,mail1,intg1,xcoor1,mode1);

% On passe en matriciel
% """""""""""""""""""""
disp('====> assemblages')

% numerotation des noeuds
numer1 = nmail1{1}.MAIL';
% liste des noms de ddl
[listDdlDualt1,listDdlPrimt1] = ListDdlModl2(modl1);
% matrices de mapping des ddl
[mapddlDualt1,mapddlPrimt1] = MapddlRigi2(modl1,mail1, ...
                                          numer1,listDdlDualt1, ...
                                          numer1,listDdlPrimt1);

%KT0 = RigiToMatrix(rigit1,modl1,mail1, ...
%                   numer1,mapddlDualt1,listDdlDualt1, ...
%                   numer1,mapddlPrimt1,listDdlPrimt1);
%ndd0 = size(KT0,1);

% On collapse les ddl bloques !
% Ddl a meme deplacement sur Y
d2 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrimt1(:,[2]),[{'BY'}]);
ddlub = sort(find(d2'));
ddlub0 = ddlub(1); % on collapse les ddl sur celui la...

U2 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrimt1,listDdlPrimt1);
ii = ddlub(2:end);
jj = setdiff([1:length(U2)],ii);
U2 = U2(jj,:);

for ii = 2:length(ddlub)
  iddl = ddlub(ii);
  [i,j] = find(mapddlPrimt1==iddl);
  for jj = 1:length(j)
    mapddlPrimt1(i(jj),j(jj)) = ddlub0;
  end
  [i,j] = find(mapddlPrimt1>iddl);
  for jj = 1:length(j)
    mapddlPrimt1(i(jj),j(jj)) = mapddlPrimt1(i(jj),j(jj)) - 1;
  end
  [i,j] = find(mapddlDualt1==iddl);
  for jj = 1:length(j)
    mapddlDualt1(i(jj),j(jj)) = ddlub0;
  end
  [i,j] = find(mapddlDualt1>iddl);
  for jj = 1:length(j)
    mapddlDualt1(i(jj),j(jj)) = mapddlDualt1(i(jj),j(jj)) - 1;
  end
  ddlub(ii+1:end) = ddlub(ii+1:end) - 1;
end

%ii = findoccur(ddlub,mapddlPrimt1(:,2)');
%mapddlPrimt1(ii,2) = ddlub0;
%ii = findoccur(ddlub,mapddlDualt1(:,2)');
%mapddlDualt1(ii,2) = ddlub0;

% Assemblage de la matrice de rigidite collapsee !
KT1 = RigiToMatrix(rigit1,modl1,mail1, ...
                   numer1,mapddlDualt1,listDdlDualt1, ...
                   numer1,mapddlPrimt1,listDdlPrimt1);
nddl = size(KT1,1);

% Ddl a deplacement impose
d1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrimt1(:,[1 2]),[{'DX'} {'DY'}]);
u1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrimt1(:,[1 2]),[{'UX1'} {'UY1'}]);
ddlud = find(d1');
ud1 = u1(ddlud,:);

% Force exterieure (resultante)
F1 = zeros(nddl,1);
F1(ddlub0,:) = 1.;


ddlui = setdiff([1:nddl],ddlud);

B_KT1 = chol(KT1(ddlui,ddlui));  % B_KT1' * B_KT1 = KT1
U1 = zeros(nddl,1);
U1(ddlui,:) = B_KT1 \ (B_KT1' \ F1(ddlui,:));
U1(ddlud,:) = 0.;

err1 = norm(U1 - U2) / norm(U2);
if (err1 > 1.e-6)
  err1
  error('PB')
end

[chpo1,nmail1] = VectToChpo2(U1,numer1,mapddlPrimt1,listDdlPrimt1);

% On sort les resultats
% """""""""""""""""""""
% Ecriture au format AVS du resultat
%file3 = 'resu_test_Assembl.avs';
%fid = fopen(file3,'w');
%  [xcoor1a,mail1a,nmail1a,chpo1a] = PrepareAVS3(xcoor1,mail1,nmail1,chpo1);
%  Lchml = [];
%  Lcara = [];
%  ListInd = lt1;
%  error1 = Write3AVS(xcoor1a,mail1a,Lchpo1,nmail1a,Lchml,Lcara,ListInd,fid)
%fclose(fid);

ListMesh2{1} = mail1;
ListChpo2{1} = chpo1;
ListnMesh2{1} = nmail1;
ListChml2{1} = [];
ListCara2{1} = [];

file3 = 'resu_test_Assembl.inp';
error1 = WriteMergeAVS(xcoor1,ListMesh2,ListChpo2,ListnMesh2, ...
                       ListChml2,ListCara2,file3);

disp('TEST PASSE AVEC SUCCES')
quit
