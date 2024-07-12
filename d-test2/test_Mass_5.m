path(path,'../matlabEF')
path(path,'../matlabEF2')
path(path,'../matlabUtils')
idim = 3;
mode1 = 'TRID';

file1 = 'test_Mass_5.inp';
nrec1 = 1; % nombre d'enregistrements
xcrit1 = -1.e-6; % critere de proximite de noeuds
[xcoor1,ListMesh1,ListChpo1,ListnMesh1,ListChml1,ListCara1,error1] = ...
  ReadMergeAVS2(xcrit1,file1,nrec1);
xcoor1 = xcoor1(:,1:idim);
mail1 = ListMesh1{1};
chpo1 = ListChpo1{1};
nmail1 = ListnMesh1{1};
% matr1 = ListChml1{1};
clear ListMesh1 ListChpo1 ListnMesh1 ListChml1 ListCara1;

% Modeles
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','MASSE',mode1,idim);
[modl2,intg2] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);

% Materiau
mater1 = ManuChml(mail1,'RHO','',1.);
mater1 = ChmlToCham(mater1,mail1,intg1);

% Matrices de masse
mass1 = Mass7(modl1,mater1,mail1,intg1,xcoor1,mode1);
% vecteur de deplacement rotation finie de solide rigide
theta1 = pi/4.3;
Q1 = [cos(theta1) -sin(theta1) 0.
      sin(theta1) cos(theta1) 0.
      0. 0. 1.];
Q2 = [cos(theta1) 0. -sin(theta1)
      sin(theta1) 0. cos(theta1)
      0. 1. 0.];
xcoor2 = xcoor1*Q1'*Q2';
mass2 = Mass7(modl1,mater1,mail1,intg1,xcoor2,mode1);
err1 = 0.;
for izo1 = 1:length(mass1)
    err1 = err1 + max(max(max(abs(mass1{izo1}.XVAL-mass2{izo1}.XVAL))));
end
err1

% Assemblage
numer1 = nmail1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
% Champ de deplacement
U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrim1,listDdlPrim1);

nbptg1 = Nptg(mail1,intg1);
Sig1 = zeros(nbptg1*6,1);
numerinv1 = InverseList(numer1,max(numer1));
Option1 = 'Green-Lagrange';
[Fi1,W1] = Fint(Sig1,mail1,xcoor1,numerinv1,mapddlPrim1, ...
                intg2,modl2, ...
                Option1,0.*U1);
xcoor2 = xcoor1 + U1(mapddlPrim1);
[Fi2,W2] = Fint(Sig1,mail1,xcoor2,numerinv1,mapddlPrim1, ...
                intg2,modl2, ...
                Option1,0.*U1);
mater2 = mater1;
ptg1 = 1;
for izo1 = 1:length(mater1)
    rho1 = mater1{izo1}{1}.XVAL;
    WW1 = W1(ptg1:ptg1+prod(size(rho1))-1,1);
    WW1 = reshape(WW1,size(rho1,2),size(rho1,1))';
    WW2 = W2(ptg1:ptg1+prod(size(rho1))-1,1);
    WW2 = reshape(WW2,size(rho1,2),size(rho1,1))';
    mater2{izo1}{1}.XVAL = (rho1 .* WW1) ./ WW2;
    ptg1 = ptg1+prod(size(rho1));
end
mass3 = Mass7(modl1,mater2,mail1,intg1,xcoor2,mode1);
err2 = 0.;
for izo1 = 1:length(mass1)
    err2 = err2 + max(max(max(abs(mass1{izo1}.XVAL-mass3{izo1}.XVAL))));
end
err2
norm(W1-W2)/norm(W1)



disp('TEST PASSE AVEC SUCCES')
quit
