% Test poutre d'Euler-Bernoulli - courbure negligee
% raccord de poutres
clear all
close all
path(path,'../matlabEF') 
path(path,'../matlabUtils') 

mode1 = 'POUT';
idim = 2;

file1 = 'test_Pout_4.inp';
nrec1 = 1; % nombre d'enregistrements
xcrit1 = -1.e-6; % critere de proximite de noeuds
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
  ReadMergeAVS2(xcrit1,file1,nrec1);
xcoor1 = xcoor1(:,1:idim);

mail1 = ListMesh1{1};
chpo1 = ListChpo1{1};
nmail1 = ListnMesh1{1};
matr1 = ListChml1{1};

clear ListMesh1 ListChpo1 ListnMesh1 ListChml1 ListCara1;

% Modele
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);

% On cree le champ par elements aux noeuds des tangentes
% mail1 est ecrase par son reoriente eventuel
xcoor2 = [100. 0.];
[chamno1,intgno1,mail1] = LocalGeo2(modl1,mail1,intg1,xcoor1,mode1,xcoor2);

% On le lisse selon un critere de brisure de tangentes
xcrit1 = 1.e6;
chamno2 = SmoothChamno(chamno1,modl1,mail1,intgno1,xcrit1);

% Materiau
matr2 = ChmlToCham(matr1,mail1,intg1);


rigi1 = Rigi9(modl1,matr2,mail1,intg1,xcoor1,mode1);
rigi2 = Rigi9(modl1,matr2,mail1,intg1,xcoor1,mode1,chamno2);
rigi3 = Rigi9(modl1,matr2,mail1,intg1,xcoor1,mode1,chamno1);


max(max(max(rigi2{1}.XVAL - rigi3{1}.XVAL)))/max(max(max(rigi2{1}.XVAL)))

% Assemblage
numer1 = RenumberRcm(mail1,xcoor1);
% nn1 = ChangeMesh2(mail1,'POI1'); numer1 = nn1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
K1 = RigiToMatrix(rigi1,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);
K2 = RigiToMatrix(rigi2,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);
K3 = RigiToMatrix(rigi3,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);
% Force imposee
F1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlDual1,listDdlDual1);
ddlf = find(F1);

% Deplacement impose
UD1 = ChpoToVect3(chpo1,nmail1,numer1,[1:length(numer1)]',[{'UD'}]);
ino1 = find(UD1);
ddld = mapddlPrim1(ino1,:);
nddl = size(F1,1);
ddli = setdiff([1:nddl],ddld);

% Resolutions
C_K1 = chol(K1(ddli,ddli));
C_K2 = chol(K2(ddli,ddli));
C_K3 = chol(K3(ddli,ddli));
U1 = zeros(nddl,1);
U1(ddli,:) = C_K1 \ (C_K1' \ F1(ddli,:)); 
U1(ddld,:) = 0.;
U2 = zeros(nddl,1);
U2(ddli,:) = C_K2 \ (C_K2' \ F1(ddli,:)); 
U2(ddld,:) = 0.;
U3 = zeros(nddl,1);
U3(ddli,:) = C_K3 \ (C_K3' \ F1(ddli,:)); 
U3(ddld,:) = 0.;

nbelt1 = Nelements2(mail1);
disp(['N1 ' int2str(nbelt1/2) ' U1 ' num2str(U1(ddlf,:)) ...
     ' U2 ' num2str(U2(ddlf,:)) ' U3 ' num2str(U3(ddlf,:))])

file2 = 'test_Pout_4U.inp';
fid = fopen(file2,'w');
  [Lchpo1,nmail1] = LvectToLchpo2([U2 U3],numer1,mapddlPrim1,listDdlPrim1);
  [xcoor1a,mail1a,nmail1a,Lchpo1a] = PrepareAVS3(xcoor1,mail1,nmail1,Lchpo1);
  Lchml = [];
  Lcara = [];
  ListInd = [1 2];
  error1 = Write3AVS(xcoor1a,mail1a,Lchpo1a,nmail1a,Lchml,Lcara,ListInd,fid)
fclose(fid);

disp('TEST PASSE AVEC SUCCES')
quit

keyboard

fid = fopen('test_Pout_4.pos','w');
  xcoor2 = [xcoor1 zeros(size(xcoor1,1),1)];
  error1 = Write1GMSH(xcoor2,mail1,chamno1,fid);
fclose(fid);
fid = fopen('test_Pout_4a.pos','w');
  xcoor2 = [xcoor1 zeros(size(xcoor1,1),1)];
  error1 = Write1GMSH(xcoor2,mail1,chamno2,fid);
fclose(fid);

% Verification
XVAL1 = U1' * K1 * U1;
if (abs(XVAL1 - 1.3349325825)/1.3349325825 > 1.E-10)
  XVAL1
  error('mauvaise valeur')
end

err1 = norm(full(K1 - K2)) / norm(full(K1));
if (err1 > 1e-7)
  err1
  error('erreur')
end

err2 = norm(full(K1 - K3)) / norm(full(K1))
err3 = norm(full(B1 - B3)) / norm(full(B1))
err4 = norm(full(bb1 - bb3)) / norm(full(bb1))
err5 = norm(full(dd1 - dd3)) / norm(full(dd1))

disp('TEST PASSE AVEC SUCCES')
quit
