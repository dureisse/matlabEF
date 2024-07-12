clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')
mode1 = 'TRID';
idim = 3;

% Recuperation des donnees
% """"""""""""""""""""""""
disp('Recuperation des donnees')
fid = fopen('test_Poreux_5.inp','rt');
  [xcoor1,mail1,chpo1,chml1,cara1,error1] = Read1AVS6(fid);
fclose(fid);
%xcoor1 = xcoor1(:,1:idim);

nmail1 = ChangeMesh2(mail1,'POI1');
numer1 = nmail1{1}.MAIL';


% Modele couple
% """""""""""""
[modl1,intg1] = ModlIntg13(mail1,'POREUX','RIGIDITE',mode1,idim);

nbelt1 = Nelements2(mail1);
nbptg1 = Nptg(mail1,intg1);
numerptg1 = [1:nbptg1];

% Coefficients materiau
% """""""""""""""""""""
mater1 = ChmlToCham(chml1,mail1,intg1);

% Matrice de rigidite/compressibilite couplee
% """""""""""""""""""""""""""""""""""""""""""
disp('Matrice de rigidite/compressibilite couplee')
[rigit1,bsigt1,bepst1] = RigiCompress8(modl1,mater1,mail1,intg1, ...
                                       xcoor1,mode1,'GenDualOp','PrimOp');

% Assemblage (sur tous les ddl : modl1)
% Vont donc ensembles :
% (numer1,mapddlDualt1,listDdlDualt1)
% (numer1,mapddlPrimt1,listDdlPrimt1)
% (numerptg1,mapcomDualt1,listComDualt1)
% (numerptg1,mapcomPrimt1,listComPrimt1)
disp('Assemblage')
[listDdlDualt1,listDdlPrimt1] = ListDdlModl2(modl1);
[listComDualt1,listComPrimt1] = ListCompModl3(modl1);
[mapddlDualt1,mapddlPrimt1] = MapddlRigi2(modl1,mail1, ...
                                          numer1,listDdlDualt1, ...
                                          numer1,listDdlPrimt1);
[mapcomDualt1,mapddlDualt0] = MapcomB2(modl1,mail1,intg1, ...
                                       numerptg1,listComDualt1, ...
                                       numer1,listDdlDualt1, ...
                                       'DUAL');
[mapcomPrimt1,mapddlPrimt0] = MapcomB2(modl1,mail1,intg1, ...
                                       numerptg1,listComPrimt1, ...
                                       numer1,listDdlPrimt1, ...
                                       'PRIMAL');
KT1 = RigiToMatrix(rigit1,modl1,mail1, ...
                   numer1,mapddlDualt1,listDdlDualt1, ...
                   numer1,mapddlPrimt1,listDdlPrimt1);
BSIGT1 = BToMatrix3(bsigt1,modl1,mail1,intg1, ...
                    numer1,mapddlDualt1,listDdlDualt1, ...
                    numerptg1,mapcomDualt1,listComDualt1, ...
                    'DUAL');
BEPST1 = BToMatrix3(bepst1,modl1,mail1,intg1, ...
                    numer1,mapddlPrimt1,listDdlPrimt1, ...
                    numerptg1,mapcomPrimt1,listComPrimt1, ...
                    'PRIMAL');


% Modele de permeabilite
% """"""""""""""""""""""
%mail2 = ChangeMesh2(mail1,'TRI3');
mail2 = mail1;
[modl2,intg2] = ModlIntg13(mail2,'FLUIDE','PERMEABILITE',mode1,idim);

% Matrice de permeabilite
% """""""""""""""""""""""
disp('Matrice de permeabilite')
[perm2,bsig2,b2] = Perm6(modl2,mater1,mail2,intg2,xcoor1,mode1, ...
                         'GenDualOp','PrimOp');

% Assemblage (sur tous les ddl : modl1)
disp('Assemblage')
HH1 = RigiToMatrix(perm2,modl2,mail2, ...
                   numer1,mapddlDualt1,listDdlDualt1, ...
                   numer1,mapddlPrimt1,listDdlPrimt1);


% CL Solides
U1 = ChpoToVect3(chpo1,nmail1,numer1,[1:length(numer1)]',[{'UD1'}]);
U2 = ChpoToVect3(chpo1,nmail1,numer1,[1:length(numer1)]',[{'UD2'}]);
U3 = ChpoToVect3(chpo1,nmail1,numer1,[1:length(numer1)]',[{'UD3'}]);
ino1 = find(U1==1);
ino2 = find(U2==1);
ino3 = find(U3==1);
% UX=UY=0
inot1 = unique([ino1' ino2' ino3']);
ddl1 = mapddlPrimt1(inot1,1:2); ddl1 = reshape(ddl1,1,prod(size(ddl1)));
% UZ=0
ddl2 = mapddlPrimt1(ino1,3); ddl2 = reshape(ddl2,1,prod(size(ddl2)));
% bloques a 0
ddl0 = unique([ddl1 ddl2]);
% UZ=impose
ddl3 = mapddlPrimt1(ino3,3); ddl3 = reshape(ddl3,1,prod(size(ddl3)));
U3 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrimt1(:,1:3),[{'UX3'} {'UY3'} {'UZ3'}]);

% ddl solides
ddlu = mapddlPrimt1(:,1:3); ddlu = reshape(ddlu,1,prod(size(ddlu)));
% ddl fluide
ddlp = mapddlPrimt1(:,4); ddlp = reshape(ddlp,1,prod(size(ddlp)));


T1 = 1.;
h = T1 / 100;
lt1 = [0.:h:T1]; theta1 = 1.;
lu1 = lt1/T1;
nddl = size(KT1,2);
ddld = [ddl0 ddl3];
ddli = setdiff([1:nddl],ddld);
ddlui = setdiff(ddlu,ddld);
Ud = U3(ddld,:);


% On suppose le pas de temps cst
M = (1./theta1/h)*KT1 - HH1;
N = KT1 + (1.-theta1)*h*HH1;
disp('====> LHS factorization')
Mii = M(ddli,ddli);
p = symrcm(Mii);
[L_M,U_M] = lu(Mii(p,p));         % Mii(p,p) = L_M * U_M
Mid = M(ddli,ddld);


X = zeros(nddl,0);
it1 = 1;
%
% Pression impose a t=0 (CI)
  Yi = zeros(nddl,1);
% Forces imposees a t=0 (CL)
  Ci = zeros(nddl,1);
  CC = Ci;
  Ci = Ci - KT1 * Yi;

  disp(['====> time step ' int2str(it1) ' resolution'])
% On calcule U(t=0)
  Ui = KT1(ddlui,ddlui) \ Ci(ddlui,:);
% On calcule X(t=0)
  Xi = Yi;
  Xi(ddlui,:) = Ui;
%
  X(:,it1) = Xi;
%
% On calcul AtX
  AtX = (1.-theta1)*h*(CC + HH1 * Xi) + KT1*Xi;

for it1 = 2:length(lt1)

  Ci = zeros(nddl,1);
  Ci(ddlu,:) = 0.;
  Ci(ddlp,:) = 0.;
  Ci = Ci + (1./theta1/h)*AtX;
  Xd = lu1(it1)*Ud; 
  Ci = Ci(ddli,:);
  Ci = Ci - Mid * Xd;

  disp(['====> time step ' int2str(it1) ' resolution'])

  Xi = zeros(nddl,1);
  sol = zeros(size(Ci));
  sol(p,:) = U_M \ (L_M \ Ci(p,:)); clear Ci;
  Xi(ddli,:) = sol; clear sol;
  Xi(ddld,:) = Xd;



  X(:,it1) = Xi;
  AtX = (1./theta1)*(KT1*Xi - (1.-theta1)*AtX);
  clear Xi;
end



[Lchpo1,nmail1] = LvectToLchpo2(X,numer1,mapddlPrimt1(:,1:3),listDdlPrimt1(1:3));
[Lchpo2,nmail2] = LvectToLchpo2(X,numer1,mapddlPrimt1(:,4),listDdlPrimt1(4));
clear Lchamno1 Lchamno2;
for i = 1:length(lt1)
  chpo1 = Lchpo1{i};
  chpo2 = Lchpo2{i};
  [chamno1,junk] = ChpoToChamno3(chpo1,nmail1,mail1);     
  [chamno2,junk] = ChpoToChamno3(chpo2,nmail2,mail1);     
  Lchamno1{i} = chamno1;
  Lchamno2{i} = chamno2;
  clear chamno1 chamno2 junk chpo1 chpo2;
end

file = 'test_Poreux_5_u.pos';
fid = fopen(file,'w');
 error1 = WriteLChGMSH3(fid,'Displacement',xcoor1,lt1, ...
                        mail1,Lchamno1,[{'UX'} {'UY'} {'UZ'}], ...
                        [],[],[])
fclose(fid)
file = 'test_Poreux_5_p.pos';
fid = fopen(file,'w');
 error1 = WriteLChGMSH3(fid,'Pore_pressure',xcoor1,lt1, ...
                        mail1,Lchamno2,[{'P'}], ...
                        [],[],[])
fclose(fid)

disp('TEST PASSE AVEC SUCCES')
quit
