% Cas TET5b

% path(path,'matlabEF')
% path(path,'matlabUtils')
clear all
close all
addpath('../matlabEF')
addpath('../matlabUtils')
mode1 = 'TRID';
idim = 3;

% Recuperation des donnees
% """"""""""""""""""""""""
disp('Recuperation des donnees')
file1 = 'test_Poreux_6.inp';
xcrit1 = 1.e-5;
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1, ...
  ListCara1,error1] = ReadMergeAVS2(xcrit1,file1,1);
mail1 = ListMesh1{1};
chpo1 = ListChpo1{1};
nmail1 = ListnMesh1{1};
chml1 = ListChml1{1};
clear ListMesh1 ListChpo1 ListnMesh1 ListChml1 ListCara1;

%fid = fopen('test_Poreux_6.avs','rt');
%  [xcoor1,mail1,chpo1,chml1,cara1,error1] = Read1AVS6(fid);
%fclose(fid);
%xcoor1 = xcoor1(:,1:idim);

% On passe en tetrahedre a 5 noeuds pour le poreux
% (noeud central de type bulle)
[mail_TET5,xcoor2] = ChangeMesh2(mail1,'TET5b',xcoor1);
xcoor1 = [xcoor1 ; xcoor2]; clear xcoor2;

nmail2 = nmail1;
nmail1 = ChangeMesh2(mail_TET5,'POI1');
numer1 = nmail1{1}.MAIL';

% Modele couple
% """""""""""""
[modl1,intg1] = ModlIntg13(mail_TET5,'POREUX','RIGIDITE',mode1,idim);

nbelt1 = Nelements2(mail_TET5);
nbptg1 = Nptg(mail_TET5,intg1);
numerptg1 = [1:nbptg1];

% Coefficients materiau
% """""""""""""""""""""
mater1 = ChmlToCham(chml1,mail_TET5,intg1);

% Matrice de rigidite/compressibilite couplee
% """""""""""""""""""""""""""""""""""""""""""
disp('Matrice de rigidite/compressibilite couplee')
% [rigit1,bsigt1,bepst1] = RigiCompress7(modl1,mater1,mail_TET5,intg1, ...
%                                        xcoor1,mode1);
[rigit1,bsigt1,bepst1] = RigiCompress8(modl1,mater1,mail_TET5,intg1, ...
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
[mapddlDualt1,mapddlPrimt1] = MapddlRigi2(modl1,mail_TET5, ...
                                          numer1,listDdlDualt1, ...
                                          numer1,listDdlPrimt1);
[mapcomDualt1,mapddlDualt0] = MapcomB2(modl1,mail_TET5,intg1, ...
                                       numerptg1,listComDualt1, ...
                                       numer1,listDdlDualt1, ...
                                       'DUAL');
[mapcomPrimt1,mapddlPrimt0] = MapcomB2(modl1,mail_TET5,intg1, ...
                                       numerptg1,listComPrimt1, ...
                                       numer1,listDdlPrimt1, ...
                                       'PRIMAL');
KT1 = RigiToMatrix(rigit1,modl1,mail_TET5, ...
                   numer1,mapddlDualt1,listDdlDualt1, ...
                   numer1,mapddlPrimt1,listDdlPrimt1);
BSIGT1 = BToMatrix3(bsigt1,modl1,mail_TET5,intg1, ...
                    numer1,mapddlDualt1,listDdlDualt1, ...
                    numerptg1,mapcomDualt1,listComDualt1, ...
                    'DUAL');
BEPST1 = BToMatrix3(bepst1,modl1,mail_TET5,intg1, ...
                    numer1,mapddlPrimt1,listDdlPrimt1, ...
                    numerptg1,mapcomPrimt1,listComPrimt1, ...
                    'PRIMAL');


% Modele de permeabilite
% """"""""""""""""""""""
mail2 = mail1;
%[modl2,intg2] = ModlIntg13(mail2,'FLUIDE','PERMEABILITE',mode1,idim);
[modl2,intg2] = ModlIntg13(mail2,'FLUIDE','PERMEABILITE',mode1,idim,intg1);

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
U1 = ChpoToVect3(chpo1,nmail2,numer1,[1:length(numer1)]',[{'UD1'}]);
U2 = ChpoToVect3(chpo1,nmail2,numer1,[1:length(numer1)]',[{'UD2'}]);
U3 = ChpoToVect3(chpo1,nmail2,numer1,[1:length(numer1)]',[{'UD3'}]);
ino1 = find(U1==1);
ino2 = find(U2==1);
ino3 = find(U3==1);
% UZ=0
ddl1 = mapddlPrimt1(ino1,3); ddl1 = reshape(ddl1,1,prod(size(ddl1)));
% UY=0
ddl2 = mapddlPrimt1(ino2,2); ddl2 = reshape(ddl2,1,prod(size(ddl2)));
% UX=0
ddl3 = mapddlPrimt1(ino3,1); ddl3 = reshape(ddl3,1,prod(size(ddl3)));
ddlud = [ddl1 ddl2 ddl3];
% Force imposee
F1 = ChpoToVect3(chpo1,nmail2,numer1,mapddlDualt1,listDdlDualt1);

% CL Fluide
P1 = ChpoToVect3(chpo1,nmail2,numer1,[1:length(numer1)]',[{'PD'}]);
ino4 = find(P1==1);
ddlpd = mapddlPrimt1(ino4,4); ddlpd = reshape(ddlpd,1,prod(size(ddlpd)));
pd1 = ChpoToVect3(chpo1,nmail2,numer1,mapddlPrimt1(:,4),[{'P1'}]);
pd1 = pd1(ddlpd,:);

% ddl solides
ddlu = mapddlPrimt1(:,1:3); ddlu = reshape(ddlu,1,prod(size(ddlu)));
% ddl fluide
ddlp = mapddlPrimt1(:,4); ddlp = reshape(ddlp,1,prod(size(ddlp)));


T1 = 1.;
T1 = 1.e-5;
h = T1 / 100;
lt1 = [0.:h:T1]; theta1 = 1.;
lf1 = lt1/T1;
lp1 = lf1;
Dlf1 = DerivateTheta2(lt1,lf1,theta1,0.);

nddl = size(KT1,2);
ddld = [ddlud ddlpd];
ddli = setdiff([1:nddl],ddld);
ddlui = setdiff(ddlu,ddlud);


% On suppose le pas de temps cst
M = (1./theta1/h)*KT1 - HH1;
%N = KT1 + (1.-theta1)*h*HH1;
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
  Yi(ddlpd,:) = lp1(it1) * pd1;
% Forces imposees a t=0 (CL)
  Ci = lf1(it1) * F1;
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
  Xpi = 0.*Xi; % derivee
  Xti = Xi + (1.-theta1)*h*Xpi;


if (1 == 0)
% On calcul AtX
  AtX = (1.-theta1)*h*CC;
  AtX = AtX(ddli,:);
  AtX = AtX + (1.-theta1)*h*HH1(ddli,ddli)*Xi(ddli,:)+KT1(ddli,ddli)*Xi(ddli,:);
%  AtX = (1.-theta1)*h*CC + (1.-theta1)*h*HH1*Xi+KT1*Xi;
end

for it1 = 2:length(lt1)

  Ci = Dlf1(it1) * F1;
%  Ci = Ci + (1./theta1/h)*AtX;
  Ci = Ci + (1./theta1/h)*KT1*Xti;
  Ci = Ci(ddli,:);

%  Ci = lf1(it1) * F1 + (1./theta1/h)*AtX;
%  Ci = Ci(ddli,:);

  Xd = zeros(length(ddld),1);
  Xd(length(ddlud)+1:length(ddld),:) = lp1(it1) * pd1;
  Ci = Ci - Mid * Xd;

  disp(['====> time step ' int2str(it1) ' resolution'])

  Xi = zeros(nddl,1);
  sol = zeros(size(Ci));
  sol(p,:) = U_M \ (L_M \ Ci(p,:)); clear Ci;
  Xi(ddli,:) = sol; clear sol;
  Xi(ddld,:) = Xd;



  X(:,it1) = Xi;

  Xti = (1./theta1)*Xi - ((1.-theta1)/theta1) * Xti;
%  AtX = (1./theta1)*(KT1(ddli,ddli)*Xi(ddli,:) - (1.-theta1)*AtX);
%  AtX = (1./theta1)*(KT1*Xi - (1.-theta1)*AtX);
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

file = 'test_Poreux_6_u.pos';
fid = fopen(file,'w');
 error1 = WriteLChGMSH3(fid,'Displacement',xcoor1,lt1, ...
                        mail1,Lchamno1,[{'UX'} {'UY'} {'UZ'}], ...
                        [],[],[])
fclose(fid)
file = 'test_Poreux_6_p.pos';
fid = fopen(file,'w');
 error1 = WriteLChGMSH3(fid,'Pore_pressure',xcoor1,lt1, ...
                        mail1,Lchamno2,[{'P'}], ...
                        [],[],[])
fclose(fid)

disp('TEST PASSE AVEC SUCCES')
quit
