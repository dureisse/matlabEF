% Test calcul forces interieures
clear all
close all

path(path,'../matlabEF') 
path(path,'../matlabEF2') 
path(path,'../matlabUtils') 

mode1 = 'TRID';
idim = 3;

file1 = 'test_Fint.inp';
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

% Modeles
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);

% Materiau
matr2 = ChmlToCham(matr1,mail1,intg1);

% Rigidite
[rigi1,bsig1] = Rigi9(modl1,matr2,mail1,intg1,xcoor1,mode1,'GenDualOp');

% Assemblage
numer1 = nmail1{1}.MAIL';
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
% K1 = RigiToMatrix(rigi1,modl1,mail1, ...
%                   numer1,mapddlDual1,listDdlDual1, ...
%                   numer1,mapddlPrim1,listDdlPrim1);

% Champ de deplacement
U1 = ChpoToVect3(chpo1,nmail1,numer1,mapddlPrim1,listDdlPrim1);

% Champ de contraintes
nbptg1 = Nptg(mail1,intg1);
numer2 = [1:nbptg1];
[listComDual1,listComPrim1] = ListCompModl3(modl1);                           
[mapcomDual1,mapddlDual0] = MapcomB2(modl1,mail1,intg1, ...
                                     numer2,listComDual1, ...
                                     numer1,listDdlDual1, ...
                                     'DUAL');
BS1 = BToMatrix3(bsig1,modl1,mail1,intg1, ...
                numer1,mapddlDual1,listDdlDual1, ...
                numer2,mapcomDual1,listComDual1, ...
                'DUAL');
nbddl = size(BS1,1);

cham1 = ChmlToCham(matr1,mail1,intg1);
Sig1 = ChamToVect3(cham1,mail1,intg1,numer2,mapcomDual1,listComDual1);

% Cas classique HPP matriciel
Fint1 = BS1' * Sig1;
Sig2 = rand(size(Sig1));
Fint2 = BS1' * Sig2;

% Sans matrice, HPP
numerinv1 = InverseList(numer1,max(numer1));
Option1 = 'HPP';
[Fi1,W1] = Fint(Sig1,mail1,xcoor1,numerinv1,mapddlPrim1, ...
           intg1,modl1, ...
           Option1,U1);
[Fi2,W2] = Fint(Sig2,mail1,xcoor1,numerinv1,mapddlPrim1, ...
           intg1,modl1, ...
           Option1,U1);
% HPP energy check Sig1*Eps1 = Fi1'*U1
Eps1 = Epsil(mail1,xcoor1,numerinv1,mapddlPrim1, ...
             intg1,modl1, ...
             Option1,U1);
V1 = U1;
Epsp1 = Epsilpoint(mail1,xcoor1,numerinv1,mapddlPrim1, ...
             intg1,modl1, ...
             Option1,U1,V1);
Err0 = max(abs(Eps1-Epsp1))/max(abs(Epsp1));
if (Err0 > 1e-10)
    error('pb hpp vitesse')
end
nbptgt1 = length(W1);
Ener1 = 0;
for ptgt1 = 1:nbptgt1
    en1 = Sig1(1+(ptgt1-1)*6:6+(ptgt1-1)*6,1)' * ...
          Epsp1(1+(ptgt1-1)*6:6+(ptgt1-1)*6,1);
    Ener1 = Ener1 + W1(ptgt1)*en1;
end
Err1 = abs(Fi1'*V1 - Ener1)/abs(Ener1);
if (Err1 > 1e-8)
    error('pb energie 1')
end
% HPP energy check Sig2*Eps1 = Fi2'*U1
Eps1 = Epsil(mail1,xcoor1,numerinv1,mapddlPrim1, ...
             intg1,modl1, ...
             Option1,U1);
V2 = U1;
nbptgt1 = length(W2);
Ener2 = 0;
for ptgt1 = 1:nbptgt1
    en1 = Sig2(1+(ptgt1-1)*6:6+(ptgt1-1)*6,1)' * ...
          Eps1(1+(ptgt1-1)*6:6+(ptgt1-1)*6,1);
    Ener2 = Ener2 + W2(ptgt1)*en1;
end
Err2 = abs(Fi2'*V2 - Ener2)/abs(Ener2);
if (Err2 > 1e-8)
    error('pb energie 2')
end
% Check
err1 = max(abs(Fint1 - Fi1)) / max(abs(Fint1));
if (err1 > 1e-8)
    error('erreur 1')
end
err1 = max(abs(Fint2 - Fi2)) / max(abs(Fint2));
if (err1 > 1e-8)
    error('erreur 2')
end


% Sans matrice, Green-Lagrange, mais petits deplacements
Sig0 = Sig2;
Sig0 = Sig1;
U2 = 1e-5*U1;
Option1 = 'Green-Lagrange';
[Fi3,W3] = Fint(Sig0,mail1,xcoor1,numerinv1,mapddlPrim1, ...
                intg1,modl1, ...
                Option1,U2);
err3 = max(abs(Fint1 - Fi3)) / max(abs(Fint1))
Fi4 = Fint(Sig0,mail1,xcoor1,numerinv1,mapddlPrim1, ...
           intg1,modl1, ...
           Option1,U1);
% Energy check Sig0*Eps4 = Fi4'*U1
Eps4 = Epsil(mail1,xcoor1,numerinv1,mapddlPrim1, ...
             intg1,modl1, ...
             Option1,U1);
V4 = U1;
Epsp4 = Epsilpoint(mail1,xcoor1,numerinv1,mapddlPrim1, ...
                   intg1,modl1, ...
                   Option1,U1,V4);
nbptgt1 = length(W3);
Ener4 = 0;
for ptgt1 = 1:nbptgt1
    en1 = Sig0(1+(ptgt1-1)*6:6+(ptgt1-1)*6,1)' * ...
          Epsp4(1+(ptgt1-1)*6:6+(ptgt1-1)*6,1);
%%%    NON : Pas Eps4... Si on prend d/dtEps4 avec V = U, ca donne un terme en plus
    Ener4 = Ener4 + W3(ptgt1)*en1;
end
Err4 = abs(Fi4'*V4 - Ener4)/abs(Ener4);
if (Err4 > 1e-8)
    Fi4'*V4
    Ener4
    error('pb energie 4')
end
% Sans matrice, Material Hencky Isotrope Elastique, mais petits deplacements
Option1 = 'Material-Hencky-IsotrElas';
% Epsp4a = Epsilpoint(mail1,xcoor1,numerinv1,mapddlPrim1, ...
%                    intg1,modl1, ...
%                    'Green-Lagrange',U2,V4);
% Epsp4b = Epsilpoint(mail1,xcoor1,numerinv1,mapddlPrim1, ...
%                    intg1,modl1, ...
%                    Option1,U2,V4);
% err3b = max(abs(Epsp4a - Epsp4b)) / max(abs(Epsp4a))
Eps4a = Epsil(mail1,xcoor1,numerinv1,mapddlPrim1, ...
              intg1,modl1, ...
              'Green-Lagrange',U2);
Eps4b = Epsil(mail1,xcoor1,numerinv1,mapddlPrim1, ...
              intg1,modl1, ...
              'Material-Hencky',U2);
Sig4b = 2*Eps4b; % ELastique isotrope, meme base propre
[Fi4a,W4a] = Fint(Sig4b,mail1,xcoor1,numerinv1,mapddlPrim1, ...
                  intg1,modl1, ...
                  'Green-Lagrange',U2);
[Fi4b,W4b] = Fint(Sig4b,mail1,xcoor1,numerinv1,mapddlPrim1, ...
                  intg1,modl1, ...
                  Option1,U2);
err30b = max(abs(Eps4a - Eps4b)) / max(abs(Eps4a))
err3b = max(abs(Fi4a - Fi4b)) / max(abs(Fi4a))

% Random fields
U5 = rand(size(U1));
V5 = rand(size(U1));
Sig5 = rand(size(Sig1));
Option1 = 'Green-Lagrange';
[Fi5,W5] = Fint(Sig5,mail1,xcoor1,numerinv1,mapddlPrim1, ...
                intg1,modl1, ...
                Option1,U5);
Epsp5 = Epsilpoint(mail1,xcoor1,numerinv1,mapddlPrim1, ...
                   intg1,modl1, ...
                   Option1,U5,V5);
nbptgt1 = length(W5);
Ener5 = 0;
for ptgt1 = 1:nbptgt1
    en1 = Sig5(1+(ptgt1-1)*6:6+(ptgt1-1)*6,1)' * ...
          Epsp5(1+(ptgt1-1)*6:6+(ptgt1-1)*6,1);
    Ener5 = Ener5 + W5(ptgt1)*en1;
end
Err5 = abs(Fi5'*V5 - Ener5)/abs(Ener5);
if (Err5 > 1e-8)
    Fi5'*V5
    Ener5
    error('pb energie 5')
end

disp('tester ausi Material-Hencky-IsotrElas')
disp('attention contrainte doit avoir meme base propre que deformation')




disp('TEST PASSE AVEC SUCCES')
quit
