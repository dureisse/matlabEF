% On teste le passage de maillage en graphe et graphe dual,
% ainsi que la numerotation locale des noeuds

clear all
close all
path(path,'../matlabEF')
path(path,'../matlabUtils')

idim = 2;

% On lit les donnees
file1 = 'test_AVS_2.inp';
xcrit1 = 1.e-6;
[xcoort1,ListMesh1,ListChpo1,ListnMesh1,ListChml1, ...
  ListCara1,error1] = ReadMergeAVS2(xcrit1,file1,2);
xcoort1 = xcoort1(:,1:idim);

mail1 = ListMesh1{2};

% Numerotation locale des noeuds
nmail1 = ChangeMesh2(mail1,'POI1');
numer1 = nmail1{1}.MAIL';
numer_inv1 = InverseList(numer1,max(numer1));
clear nmail1;
mail2 = RenumMesh(mail1,numer_inv1);

% On tranforme mail2 en graphe
[lnod1,lind1] = MeshToGraph(mail2);
% On dualise
[lnod2,lind2] = InverseGraph2(lnod1,lind1);
% On dualise
[lnod3,lind3] = InverseGraph2(lnod2,lind2);
% L'ordre des noeuds n'est pas forcement respecte

err1 = 0;
for el1 = 1:length(lind1)-1
  err1 = norm(sort(lnod1(lind1(el1):lind1(el1+1)-1)) - ...
              sort(lnod3(lind3(el1):lind3(el1+1)-1)));
end
err1 = err1/norm(lnod1);
err2 = norm(lind1 - lind3)/norm(lind1);
if (err1 > 1.e-10) | (err2 > 1.e-10)
  err1
  err2
  error('PB')
end

% On dualise
[lnod2,lind2,lord2] = InverseGraph2(lnod1,lind1);
% On dualise
[lnod3,lind3] = InverseGraph2(lnod2,lind2,lord2);
% L'ordre des noeuds doit etre respecte
err1a = norm(lnod1-lnod3)/norm(lnod1);
err2a = norm(lind1-lind3)/norm(lind1);
if (err1a > 1.e-10) | (err2a > 1.e-10)
  err1a
  err2a
  error('PB bis')
end

% Retour a la numerotation globale
mail1 = RenumMesh(mail2,numer1);
ListMesh1{2} = mail1;

% On ecrit les resultats
file2 = 'test_AVS_2b.inp';
error1 = WriteMergeAVS(xcoort1,ListMesh1,ListChpo1,ListnMesh1, ...
                       ListChml1,ListCara1,file2);


% Verifications a la main
% 4---3
% |   | \
% 5---1--2
lnod1 = [2 1 3 1 5 4 3];
lind1 = [1 4 8];
[lnod2,lind2] = InverseGraph2(lnod1,lind1);
[lnod3,lind3] = InverseGraph2(lnod2,lind2);
% (lnod3,lind3) n'a pas le meme ordre que (lno1,lind1)
[lnod2,lind2,lord2] = InverseGraph2(lnod1,lind1);
[lnod4,lind4] = InverseGraph2(lnod2,lind2,lord2);
% (lnod4,lind4) a le meme ordre que (lno1,lind1)
[lnod5,lind5] = InverseGraph2(lnod4,lind4);

err4 = norm(lnod4-lnod1)/norm(lnod1) + norm(lind4-lind1)/norm(lind1);
if (err4 > 1.e-10)
  err4
  error('PB 2')
end
err5 = norm(lnod5-lnod2)/norm(lnod2) + norm(lind5-lind2)/norm(lind2);
if (err5 > 1.e-10)
  err5
  error('PB 3')
end




disp('FAIRE UN diff test_AVS_2.inp test_AVS_2b.inp')
system('diff -b test_AVS_2.inp test_AVS_2b.inp')
quit
