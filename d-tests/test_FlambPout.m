% test de flambage de poutre sur appui elastique reparti
% Le quotient de Rayleigh est alors 
% R = (0.5 * U' * K * U + 0.5 * U' * M1 * U) / (0.5 * U' * M2 * U)
% R = 0.5 * U' * (K1 + M1) * U / 0.5 * U' * M2 * U
% dans M1, rho.S = k et rho.I = 0
% dans M2, rho.S = 0 et rho.I = 1
% Ce qui donne le pb aux valeurs propres
% (K1 + M1) U - R M2 U = 0
% R est la charge critique, U le mode de flambement
%
clear all
close all
path(path,'../matlabEF') 
path(path,'../matlabUtils') 

mode1 = 'POUT';
idim = 2;

file1 = 'test_FlambPout.inp';
nrec1 = 1; % nombre d'enregistrements
xcrit1 = -1.e-6; % critere de proximite de noeuds
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
  ReadMergeAVS2(xcrit1,file1,nrec1);
xcoor1 = xcoor1(:,1:idim);

mail1 = ListMesh1{1};
nmail1 = ChangeMesh2(mail1,'POI1');
mate1 = ListChml1{1};

% Modele
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);
[modl2,intg2] = ModlIntg13(mail1,'ELASTIQUE','MASSE',mode1,idim);

% Materiau
matr1 = ChmlToCham(mate1,mail1,intg1);
mate1a = mate1; mate1a{4}.XVAL = 0. * mate1a{4}.XVAL;
mate1b = mate1; mate1b{3}.XVAL = 0. * mate1a{3}.XVAL;
matr2a = ChmlToCham(mate1a,mail1,intg2);
matr2b = ChmlToCham(mate1b,mail1,intg2);

rigi1 = Rigi9(modl1,matr1,mail1,intg1,xcoor1,mode1);
mass1 = Mass7(modl2,matr2a,mail1,intg2,xcoor1,mode1);
mass2 = Mass7(modl2,matr2b,mail1,intg2,xcoor1,mode1);

% Assemblage
nmail1 = ChangeMesh2(mail1,'POI1');
numer1 = nmail1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[listComDual1,listComPrim1] = ListCompModl3(modl1);
[listComDual2,listComPrim2] = ListCompModl3(modl2);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
K1 = RigiToMatrix(rigi1,modl1,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);
M1 = RigiToMatrix(mass1,modl2,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);
M2 = RigiToMatrix(mass2,modl2,mail1, ...
                  numer1,mapddlDual1,listDdlDual1, ...
                  numer1,mapddlPrim1,listDdlPrim1);
ddlu = mapddlPrim1(:,1)';
nddl = size(K1,2);
ddlv = setdiff([1:nddl],ddlu);
KK = K1(ddlv,ddlv) + M1(ddlv,ddlv);
MM = M2(ddlv,ddlv);
%+ 1.E-10*M1(ddlv,ddlv);
MM = 0.5 * (MM + MM');
[V,D,FLAG] = eigs(KK,MM,3,'SM');
if FLAG
  error('PB')
end
disp(['Forces critiques ' num2str(D(1,1)) ' N ' num2str(D(2,2)) ' N '])
if (1 == 0)
ll = [1:size(V,1)/2];
plot(ll,V([1:2:end],1)','b-',ll,V([1:2:end],2)','r--')
legend('mode 1','mode 2')
end

quit



