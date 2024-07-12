% Test calcul forces interieures
% grand maillage : librairie git LibHip https://github.com/diku-dk/libhip
clear all
close all

path(path,'../matlabEF') 
path(path,'../matlabEF2') 
path(path,'../matlabUtils') 

mode1 = 'TRID';
idim = 3;

m1_wg_legL

xcoor1 = msh.POS;
nbno = size(xcoor1,1);
nbddl = idim*nbno;
clear maile1; maile1 = struct('TYPE','TET4','MAIL',msh.TETS(:,1:4));
clear mail1; mail1{1} = maile1;
clear msh;

% Modeles
[modl1,intg1] = ModlIntg13(mail1,'ELASTIQUE','RIGIDITE',mode1,idim);

% Assemblage
numer1 = [1:nbno];
[listDdlDual1,listDdlPrim1] = ListDdlModl2(modl1);
[mapddlDual1,mapddlPrim1] = MapddlRigi2(modl1,mail1, ...
                                        numer1,listDdlDual1, ...
                                        numer1,listDdlPrim1);
% Champ de deplacement
U1 = rand(nbddl,1);

% Champ de contraintes
nbptg1 = Nptg(mail1,intg1);
numer2 = [1:nbptg1];
[listComDual1,listComPrim1] = ListCompModl3(modl1);                           
[mapcomDual1,mapddlDual0] = MapcomB2(modl1,mail1,intg1, ...
                                     numer2,listComDual1, ...
                                     numer1,listDdlDual1, ...
                                     'DUAL');
nbcomp = nbptg1*6;
Sig1 = rand(nbcomp,1);

% Sans matrice, HPP
numerinv1 = InverseList(numer1,max(numer1));
Option1 = 'HPP';
tic
[Fi1,W1] = Fint(Sig1,mail1,xcoor1,numerinv1,mapddlPrim1, ...
           intg1,modl1, ...
           Option1,U1);
toc
% HPP energy check Sig1*Eps1 = Fi1'*U1
Eps1 = Epsil(mail1,xcoor1,numerinv1,mapddlPrim1, ...
             intg1,modl1, ...
             Option1,U1);
V1 = U1;
Epsp1 = Epsilpoint(mail1,xcoor1,numerinv1,mapddlPrim1, ...
                   intg1,modl1, ...
                   Option1,U1,V1);
Err0 = max(abs(Epsp1-Epsp1))/max(abs(Epsp1));
if (Err0 > 1e-10)
    error('pb hpp vitesse')
end
nbptgt1 = length(W1);
Ener1 = 0;
for ptgt1 = 1:nbptgt1
    en1 = Sig1(1+(ptgt1-1)*6:6+(ptgt1-1)*6,1)' * ...
          Eps1(1+(ptgt1-1)*6:6+(ptgt1-1)*6,1);
    Ener1 = Ener1 + W1(ptgt1)*en1;
end
Err1 = abs(Fi1'*V1 - Ener1)/abs(Ener1);
if (Err1 > 1e-8)
    error('pb energie 1')
end

% Sans matrice, Green-Lagrange
Option1 = 'Green-Lagrange';
tic
[Fi5,W5] = Fint(Sig1,mail1,xcoor1,numerinv1,mapddlPrim1, ...
                intg1,modl1, ...
                Option1,U1);
toc
tic
Epsp5 = Epsilpoint(mail1,xcoor1,numerinv1,mapddlPrim1, ...
                   intg1,modl1, ...
                   Option1,U1,V1);
toc
nbptgt1 = length(W5);
Ener5 = 0;
for ptgt1 = 1:nbptgt1
    en1 = Sig1(1+(ptgt1-1)*6:6+(ptgt1-1)*6,1)' * ...
          Epsp5(1+(ptgt1-1)*6:6+(ptgt1-1)*6,1);
    Ener5 = Ener5 + W5(ptgt1)*en1;
end
Err5 = abs(Fi5'*V1 - Ener5)/abs(Ener5);
if (Err5 > 1e-8)
    Fi5'*V5
    Ener5
    error('pb energie 5')
end






disp('TEST PASSE AVEC SUCCES')
quit
